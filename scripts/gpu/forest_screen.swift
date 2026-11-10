import Foundation
import Metal

func clockTime()->Double {Double(DispatchTime.now().uptimeNanoseconds)*1e-9}
func require(_ value:Bool,_ message:String)throws{
    if !value {throw NSError(domain:"forest-screen",code:1,userInfo:[NSLocalizedDescriptionKey:message])}
}
func run()throws{
    let args=CommandLine.arguments
    try require(args.count==3,"usage: forest_screen TRACE SHADER")
    let data=try Data(contentsOf:URL(fileURLWithPath:args[1]))
    try require(data.count>=8,"truncated trace")
    let d=Int(data.withUnsafeBytes{$0.loadUnaligned(fromByteOffset:0,as:UInt32.self).littleEndian})
    let count=Int(data.withUnsafeBytes{$0.loadUnaligned(fromByteOffset:4,as:UInt32.self).littleEndian})
    let stride=6*d+5
    try require((9...14).contains(d)&&count>0&&data.count==8+stride*count,"trace bounds")
    guard let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue(),
        let input=device.makeBuffer(length:stride*count,options:.storageModeShared),
        let output=device.makeBuffer(length:count,options:.storageModeShared) else {
        throw NSError(domain:"Metal unavailable",code:1)
    }
    data.withUnsafeBytes{memcpy(input.contents(),$0.baseAddress!.advanced(by:8),stride*count)}
    let lib=try device.makeLibrary(source:String(contentsOfFile:args[2],encoding:.utf8),options:nil)
    guard let function=lib.makeFunction(name:"forest_screen")else{throw NSError(domain:"shader",code:1)}
    let pipeline=try device.makeComputePipelineState(function:function)
    let rows=input.contents().assumingMemoryBound(to:UInt8.self)
    var cpuOut=[UInt8](repeating:0,count:count)
    func cpu(){
        for i in 0..<count{
            let row=rows.advanced(by:i*stride);let c=row.advanced(by:5*d)
            let ia=Int(row[6*d]),ib=Int(row[6*d+1]),ic=Int(row[6*d+2]),id=Int(row[6*d+3])
            var result:UInt8=0
            for j in 0..<5{
                let left=j<2 ? ia+ic : ib+ic,right=j<2 ? ib+id : ia+id
                if left>=d && right>=d{result=UInt8(j+1);break}
                let a=row.advanced(by:j*d);var seen:UInt64=0;var cycles=0
                for root in 0..<d{
                    if seen&(UInt64(1)<<root) != 0{continue}
                    cycles+=1;var x=root
                    while seen&(UInt64(1)<<x)==0{seen|=UInt64(1)<<x;x=Int(a[Int(c[x])])}
                }
                if d-cycles != left && d-cycles != right{result=UInt8(j+1);break}
            }
            cpuOut[i]=result
        }
    }
    func gpu()throws{
        guard let command=queue.makeCommandBuffer(),let enc=command.makeComputeCommandEncoder()else{throw NSError(domain:"encoder",code:1)}
        enc.setComputePipelineState(pipeline);enc.setBuffer(input,offset:0,index:0);enc.setBuffer(output,offset:0,index:1)
        var params=[UInt32(d),UInt32(count)]
        params.withUnsafeMutableBytes{enc.setBytes($0.baseAddress!,length:$0.count,index:2)}
        enc.dispatchThreads(MTLSize(width:count,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:min(128,pipeline.maxTotalThreadsPerThreadgroup),height:1,depth:1))
        enc.endEncoding();command.commit();command.waitUntilCompleted()
        try require(command.status == .completed,"Metal dispatch failed")
    }
    cpu();try gpu()
    let gpuRows=output.contents().assumingMemoryBound(to:UInt8.self)
    for i in 0..<count{try require(cpuOut[i]==rows[i*stride+6*d+4]&&gpuRows[i]==cpuOut[i],"literal Rust/Swift/Metal mismatch at \(i)")}
    print("stage\tdegree\trepeat\tvariant\tseconds\tchecksum\tunits")
    for repeatID in 1...7{
        let variants=repeatID%2==0 ? ["cpu","gpu","gpu_shadow"] : ["gpu_shadow","gpu","cpu"]
        for variant in variants{
            let start=clockTime()
            if variant=="cpu"{cpu()}
            else{
                try gpu()
                if variant=="gpu_shadow"{
                    cpu()
                    try require(cpuOut.withUnsafeBytes{memcmp($0.baseAddress!,output.contents(),count)}==0,"shadow mismatch")
                }
            }
            let elapsed=clockTime()-start
            print("metal_screen\t\(d)\t\(repeatID)\t\(variant)\t\(elapsed)\t\(count)\t\(count)")
        }
    }
    print("# full_batch_literal_verified=true device=\(device.name)")
}
do{try run()}catch{fputs("forest Metal: \(error)\n",stderr);exit(1)}
