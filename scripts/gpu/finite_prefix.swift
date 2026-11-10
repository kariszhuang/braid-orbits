import Foundation
import Metal
import Darwin
typealias CPU = @convention(c)(UInt32,UInt32,UInt32,UnsafeRawPointer,UnsafeRawPointer,UnsafeMutableRawPointer)->Int32
func now()->Double{Double(DispatchTime.now().uptimeNanoseconds)*1e-9}
func require(_ b:Bool,_ s:String)throws{if !b{throw NSError(domain:"finite-prefix",code:1,userInfo:[NSLocalizedDescriptionKey:s])}}
func run()throws{
    let a=CommandLine.arguments
    try require(a.count==8,"usage: finite-prefix LIB SHADER CFILE FFILE JSON REPEATS WORKERS")
    let repeats=Int(a[6])!,workers=UInt32(a[7])!
    try require(repeats>=3&&repeats<=11&&workers>=1&&workers<=4,"bounded benchmark budget")
    guard let dev=MTLCreateSystemDefaultDevice(),let queue=dev.makeCommandQueue(),let lib=dlopen(a[1],RTLD_NOW|RTLD_LOCAL),let sym=dlsym(lib,"finite_prefix")else{throw NSError(domain:"runtime",code:1)}
    defer{dlclose(lib)}
    let cpu=unsafeBitCast(sym,to:CPU.self),begin=now()
    let source=try String(contentsOfFile:a[2],encoding:.utf8)
    let shaders=try dev.makeLibrary(source:source,options:nil)
    guard let function=shaders.makeFunction(name:"prefix")else{throw NSError(domain:"function",code:1)}
    let pipeline=try dev.makeComputePipelineState(function:function)
    let cs=try Data(contentsOf:URL(fileURLWithPath:a[3])),fs=try Data(contentsOf:URL(fileURLWithPath:a[4]))
    let count=cs.count/8,frames=fs.count/40,total=count*frames
    try require(cs.count==count*8&&fs.count==frames*40&&total>0&&total<=Int(UInt32.max),"input sizes")
    try require(2*total<=dev.maxBufferLength,"output exceeds Metal buffer limit")
    guard let cb=cs.withUnsafeBytes({dev.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)}),
          let fb=fs.withUnsafeBytes({dev.makeBuffer(bytes:$0.baseAddress!,length:$0.count,options:.storageModeShared)}),
          let co=dev.makeBuffer(length:2*total,options:.storageModeShared),let go=dev.makeBuffer(length:2*total,options:.storageModeShared)else{throw NSError(domain:"buffers",code:1)}
    let setup=now()-begin
    func cpuCall()throws->Double{let t=now();try require(cpu(UInt32(frames),UInt32(count),workers,fb.contents(),cb.contents(),co.contents())==0,"CPU failure");return now()-t}
    func gpuCall(_ group:Int)throws->(Double,Double){let t=now()
        guard let command=queue.makeCommandBuffer(),let enc=command.makeComputeCommandEncoder()else{throw NSError(domain:"encoder",code:1)}
        enc.setComputePipelineState(pipeline);enc.setBuffer(cb,offset:0,index:0);enc.setBuffer(fb,offset:0,index:1);enc.setBuffer(go,offset:0,index:2)
        var params=[UInt32(count),UInt32(total)];params.withUnsafeMutableBytes{enc.setBytes($0.baseAddress!,length:$0.count,index:3)}
        enc.dispatchThreads(MTLSize(width:total,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:min(group,pipeline.maxTotalThreadsPerThreadgroup),height:1,depth:1))
        enc.endEncoding();command.commit();command.waitUntilCompleted()
        try require(command.status == .completed,"Metal execution failure")
        return(now()-t,command.gpuEndTime-command.gpuStartTime)
    }
    var rows=[[String:Any]]()
    for group in [128,256]{_=try cpuCall();_=try gpuCall(group);try require(memcmp(co.contents(),go.contents(),2*total)==0,"warm CPU/GPU mismatch")
        for repeatIndex in 0..<repeats{var c=0.0,g=0.0,d=0.0
            for variant in (repeatIndex%2==0 ? [0,1]:[1,0]){if variant==0{c=try cpuCall()}else{let result=try gpuCall(group);g=result.0;d=result.1}}
            try require(memcmp(co.contents(),go.contents(),2*total)==0,"CPU/GPU integer mismatch")
            rows.append(["group":group,"repeat":repeatIndex,"cpu_seconds":c,"gpu_wall_seconds":g,"gpu_device_seconds":d])
        }
    }
    let result:[String:Any]=["device":dev.name,"frames":frames,"completions":count,"states":total,"cpu_workers":workers,"setup_seconds":setup,"all_integer_outputs_equal":true,"scope":"resident five-word screen; NOT full search or independent certification", "measurements":rows]
    let data=try JSONSerialization.data(withJSONObject:result,options:[.prettyPrinted,.sortedKeys])
    try data.write(to:URL(fileURLWithPath:a[5]),options:.withoutOverwriting)
    print(String(data:data,encoding:.utf8)!)
}
do{try run()}catch{fputs("\(error)\n",stderr);exit(1)}
