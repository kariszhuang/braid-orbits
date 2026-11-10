import Foundation
import Metal

func clockTime()->Double {Double(DispatchTime.now().uptimeNanoseconds)*1e-9}
func require(_ value:Bool,_ message:String)throws {
    if !value {throw NSError(domain:"middle-screen",code:1,userInfo:[NSLocalizedDescriptionKey:message])}
}
func run()throws {
    let args=CommandLine.arguments
    try require(args.count==3,"usage: forest_middle_screen TRACE SHADER")
    let data=try Data(contentsOf:URL(fileURLWithPath:args[1]))
    try require(data.count>=8,"truncated trace")
    let d=Int(data.withUnsafeBytes{$0.loadUnaligned(fromByteOffset:0,as:UInt32.self).littleEndian})
    let count=Int(data.withUnsafeBytes{$0.loadUnaligned(fromByteOffset:4,as:UInt32.self).littleEndian})
    let stride=5*d+5
    try require((9...14).contains(d)&&count>0&&data.count==8+stride*count,"trace bounds")
    guard let device=MTLCreateSystemDefaultDevice(),let queue=device.makeCommandQueue(),
        let input=device.makeBuffer(length:stride*count,options:.storageModeShared),
        let output=device.makeBuffer(length:count,options:.storageModeShared)else{throw NSError(domain:"Metal unavailable",code:1)}
    let lib=try device.makeLibrary(source:String(contentsOfFile:args[2],encoding:.utf8),options:nil)
    guard let fn=lib.makeFunction(name:"forest_middle_screen")else{throw NSError(domain:"shader",code:1)}
    let pipeline=try device.makeComputePipelineState(function:fn)
    let rows=input.contents().assumingMemoryBound(to:UInt8.self)
    let gpuRows=output.contents().assumingMemoryBound(to:UInt8.self)
    var cpuOut=[UInt8](repeating:0,count:count)
    func copyInput(){_ = data.withUnsafeBytes{memcpy(input.contents(),$0.baseAddress!.advanced(by:8),stride*count)}}
    func matches(_ row:UnsafePointer<UInt8>,_ forward:Bool)->Bool {
        let a=row,b=row.advanced(by:d),bi=row.advanced(by:2*d),c=row.advanced(by:3*d),ci=row.advanced(by:4*d)
        let expected=d-Int(row[5*d])-Int(row[5*d+1]);var cycles=0;var remaining=(UInt32(1)<<d)-1
        while remaining != 0 {
            cycles+=1;if cycles>expected{return false}
            var x=remaining.trailingZeroBitCount
            while remaining&(UInt32(1)<<x) != 0 {
                remaining &= ~(UInt32(1)<<x)
                let y=forward ? b[Int(c[Int(b[Int(ci[Int(bi[x])])])])]:ci[Int(b[Int(c[x])])]
                x=Int(a[Int(y)])
            }
            if remaining.nonzeroBitCount<expected-cycles{return false}
        }
        return cycles==expected
    }
    func cpu(){for i in 0..<count {
        let row=UnsafePointer(rows.advanced(by:i*stride))
        cpuOut[i] = !matches(row,true) ? 1 : !matches(row,false) ? 2 : 0
    }}
    func gpu()throws {
        guard let cmd=queue.makeCommandBuffer(),let enc=cmd.makeComputeCommandEncoder()else{throw NSError(domain:"encoder",code:1)}
        enc.setComputePipelineState(pipeline);enc.setBuffer(input,offset:0,index:0);enc.setBuffer(output,offset:0,index:1)
        var params=[UInt32(d),UInt32(count)]
        params.withUnsafeMutableBytes{enc.setBytes($0.baseAddress!,length:$0.count,index:2)}
        enc.dispatchThreads(MTLSize(width:count,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:min(128,pipeline.maxTotalThreadsPerThreadgroup),height:1,depth:1))
        enc.endEncoding();cmd.commit();cmd.waitUntilCompleted();try require(cmd.status == .completed,"Metal dispatch failed")
    }
    func verifyNegativeWitnesses()throws {
        for i in 0..<count {
            let code=gpuRows[i]
            try require(code<=2,"invalid GPU witness code")
            if code==0{continue} // Kept for the CPU BFS; never a positive proof.
            try require(!matches(UnsafePointer(rows.advanced(by:i*stride)),code==1),"false GPU rejection at \(i)")
        }
    }
    copyInput();cpu();try gpu();try verifyNegativeWitnesses()
    for i in 0..<count{try require(cpuOut[i]==rows[i*stride+5*d+4]&&gpuRows[i]==cpuOut[i],"Rust/Swift/Metal mismatch at \(i)")}
    print("stage\tdegree\trepeat\tvariant\tseconds\tchecksum\tunits")
    for repeatID in 1...7 {
        for variant in repeatID%2==0 ? ["cpu","gpu","gpu_checked"]:["gpu_checked","gpu","cpu"] {
            let start=clockTime()
            if variant=="cpu"{cpu()}else{copyInput();try gpu();if variant=="gpu_checked"{try verifyNegativeWitnesses()}}
            print("middle_square\t\(d)\t\(repeatID)\t\(variant)\t\(clockTime()-start)\t\(count)\t\(count)")
        }
    }
    print("# full_batch_differential=true every_gpu_rejection_cpu_verified=true device=\(device.name)")
}
do{try run()}catch{fputs("Middle Metal: \(error)\n",stderr);exit(1)}
