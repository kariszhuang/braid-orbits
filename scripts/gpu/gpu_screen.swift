import Foundation
import Metal
import Darwin

typealias Generate = @convention(c) (UInt32,UInt32,UInt32,UInt64,UnsafeMutableRawPointer)->UInt32
typealias Screen = @convention(c) (UInt32,UInt32,UInt32,UnsafeRawPointer,UnsafeMutableRawPointer)->UInt32
typealias Compact = @convention(c) (UInt32,UnsafeRawPointer,UnsafeMutableRawPointer)->UInt32
func now()->Double {Double(DispatchTime.now().uptimeNanoseconds)*1e-9}
func require(_ condition:Bool,_ message:String) throws {
    if !condition {throw NSError(domain:"gpu-screen",code:1,userInfo:[NSLocalizedDescriptionKey:message])}
}

func run() throws {
    let args=CommandLine.arguments
    try require(args.count==7,"usage: gpu_screen LIBRARY SHADER OUT ORACLE_INPUT ORACLE_OUTPUT REPEATS")
    let repeats=Int(args[6])!;try require(repeats>=3 && repeats<=31,"require 3..31 repeats")
    guard let handle=dlopen(args[1],RTLD_NOW|RTLD_LOCAL) else {throw NSError(domain:"dlopen",code:1)}
    defer {dlclose(handle)}
    guard let genSymbol=dlsym(handle,"bench_generate"),let screenSymbol=dlsym(handle,"bench_screen"),let compactSymbol=dlsym(handle,"bench_compact") else {throw NSError(domain:"dlsym",code:1)}
    let generate=unsafeBitCast(genSymbol,to:Generate.self)
    let cpuScreen=unsafeBitCast(screenSymbol,to:Screen.self)
    let compact=unsafeBitCast(compactSymbol,to:Compact.self)
    guard let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue() else {throw NSError(domain:"Metal",code:1)}
    let library=try device.makeLibrary(source:String(contentsOfFile:args[2],encoding:.utf8),options:nil)
    guard let function=library.makeFunction(name:"screen") else {throw NSError(domain:"shader",code:1)}
    let pipeline=try device.makeComputePipelineState(function:function)
    let maxN=262144,maxD=34
    guard let input=device.makeBuffer(length:maxN*3*maxD,options:.storageModeShared),
          let gpuOut=device.makeBuffer(length:maxN*4,options:.storageModeShared),
          let cpuOut=device.makeBuffer(length:maxN*4,options:.storageModeShared),
          let compactOut=device.makeBuffer(length:maxN*4,options:.storageModeShared) else {throw NSError(domain:"buffers",code:1)}

    func gpu(_ d:Int,_ n:Int,_ mode:Int) throws -> (Double,Double) {
        let start=now()
        guard let command=queue.makeCommandBuffer(),let encoder=command.makeComputeCommandEncoder() else {throw NSError(domain:"encoder",code:1)}
        encoder.setComputePipelineState(pipeline);encoder.setBuffer(input,offset:0,index:0);encoder.setBuffer(gpuOut,offset:0,index:1)
        let params=[UInt32(d),UInt32(n),UInt32(mode)]
        params.withUnsafeBytes {encoder.setBytes($0.baseAddress!,length:$0.count,index:2)}
        encoder.dispatchThreads(MTLSize(width:n,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:min(128,pipeline.maxTotalThreadsPerThreadgroup),height:1,depth:1))
        encoder.endEncoding();command.commit();command.waitUntilCompleted()
        try require(command.status == .completed,"GPU failure: \(String(describing:command.error))")
        return (now()-start,command.gpuEndTime-command.gpuStartTime)
    }
    func equal(_ n:Int) throws {try require(memcmp(cpuOut.contents(),gpuOut.contents(),n*4)==0,"CPU/GPU integer result mismatch")}
    func cpu(_ d:Int,_ n:Int,_ mode:Int) throws ->Double {
        let start=now();try require(cpuScreen(UInt32(d),UInt32(n),UInt32(mode),input.contents(),cpuOut.contents())==0,"CPU screen failed")
        return now()-start
    }
    func fill(_ d:Int,_ n:Int,_ family:Int,_ seed:UInt64) throws {
        try require(generate(UInt32(d),UInt32(n),UInt32(family),seed,input.contents())==0,"Generator failed")
    }
    // Independent Python oracle: every ordered triple of S4 permutations.
    let oracleIn=try Data(contentsOf:URL(fileURLWithPath:args[4]))
    let oracleOut=try Data(contentsOf:URL(fileURLWithPath:args[5]))
    let oracleN=oracleIn.count/12
    try require(oracleIn.count==oracleN*12 && oracleN<=maxN && oracleOut.count==oracleN*4,"invalid oracle sizes")
    oracleIn.withUnsafeBytes {memcpy(input.contents(),$0.baseAddress!,$0.count)}
    _=try cpu(4,oracleN,0);_=try gpu(4,oracleN,0);try equal(oracleN)
    try oracleOut.withUnsafeBytes {try require(memcmp(cpuOut.contents(),$0.baseAddress!,$0.count)==0,"independent exhaustive S4 oracle mismatch")}
    _=try cpu(4,oracleN,1);_=try gpu(4,oracleN,1);try equal(oracleN)
    for d in [18,25] {
        try fill(d,16,2,20260929);_=try cpu(d,16,1);_=try gpu(d,16,1);try equal(16)
        try require(cpuOut.contents().assumingMemoryBound(to:UInt32.self)[0]==1,"known weak certificate rejected")
    }
    print("GPU: \(device.name); exhaustive S4 oracle and known weak controls passed")
    var rows=["degree\tfamily\tbatch\tmode\trepeat\tcpu_seconds\tgpu_wall_seconds\tgpu_device_seconds\tcpu_pipeline_seconds\tgpu_pipeline_seconds\tsurvivors\tproducer_seconds"]
    for d in [9,10,11,18,25,34] {
        for family in (d==18 || d==25 ? [0,1,2] : [0,1]) {
            for n in [256,4096,65536,262144] {
                try fill(d,n,family,20260929)
                for mode in [0,1] {
                    // Warm resident paths; compilation and allocation excluded.
                    _=try cpu(d,n,mode);_=try gpu(d,n,mode);try equal(n)
                    for rep in 0..<repeats {
                        var c=0.0,g=0.0,deviceTime=0.0,cp=0.0,gp=0.0,producer=0.0
                        let variants=rep%2==0 ? [0,1] : [1,0]
                        for variant in variants {
                            if variant==0 {c=try cpu(d,n,mode)}
                            else {let times=try gpu(d,n,mode);g=times.0;deviceTime=times.1}
                        }
                        try equal(n)
                        var survivors:UInt32=0
                        if mode==1 {
                            let seed=UInt64(20260929+rep*3571)
                            var counts=[UInt32]()
                            for variant in variants {
                                let start=now(),generationStart=start
                                try fill(d,n,family,seed);let generated=now()-generationStart
                                if variant==0 {_=try cpu(d,n,mode)} else {_=try gpu(d,n,mode)}
                                let count=compact(UInt32(n),variant==0 ? cpuOut.contents() : gpuOut.contents(),compactOut.contents())
                                let elapsed=now()-start;counts.append(count)
                                if variant==0 {cp=elapsed;producer=generated} else {gp=elapsed}
                            }
                            try equal(n);try require(counts[0]==counts[1],"survivor count mismatch");survivors=counts[0]
                        }
                        rows.append("\(d)\t\(family)\t\(n)\t\(mode)\t\(rep)\t\(c)\t\(g)\t\(deviceTime)\t\(cp)\t\(gp)\t\(survivors)\t\(producer)")
                    }
                }
                print("passed degree=\(d) family=\(family) batch=\(n)")
                try (rows.joined(separator:"\n")+"\n").write(toFile:args[3],atomically:true,encoding:.utf8)
            }
        }
    }
}
do {try run()} catch {fputs("\(error)\n",stderr);exit(1)}
