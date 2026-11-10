// Shared discovery/benchmark pipeline; bounded shared memory, exact CPU checks.
import Foundation
import Metal
import Darwin
import CryptoKit

typealias OpenSession = @convention(c) (UnsafePointer<CChar>,UnsafePointer<CChar>)->UnsafeMutableRawPointer?
typealias CloseSession = @convention(c) (UnsafeMutableRawPointer)->Void
typealias Fill = @convention(c) (UnsafeRawPointer,UInt32,UInt32,UInt32,UInt64,UInt32,UnsafeMutableRawPointer)->UInt32
typealias ScreenCPU = @convention(c) (UInt32,UInt32,UInt32,UnsafeRawPointer,UnsafeMutableRawPointer)->UInt32
typealias Validate = @convention(c) (UnsafeRawPointer,UInt32,UInt32,UInt32,UInt32,UInt64,UnsafeRawPointer,UnsafeRawPointer,UnsafeMutableRawPointer)->UInt32
typealias Diagnostic = @convention(c) (UInt32,UInt32,UInt32,UnsafeRawPointer,UnsafeMutableRawPointer)->UInt32
func now()->Double {Double(DispatchTime.now().uptimeNanoseconds)*1e-9}
func check(_ ok:Bool,_ text:String) throws {if !ok {throw NSError(domain:"hybrid-pilot",code:1,userInfo:[NSLocalizedDescriptionKey:text])}}
func digest(_ bytes:Data)->String {SHA256.hash(data:bytes).map {String(format:"%02x",$0)}.joined()}
func packedFlags(_ flags:UnsafePointer<UInt32>,_ n:Int)throws->Data {
    var packed=Data(repeating:0,count:(n+7)/8)
    for i in 0..<n {try check(flags[i]<=1,"non-binary screening flag");if flags[i]==1 {packed[i/8]|=UInt8(1<<(i%8))}}
    return packed
}

func run() throws {
    let a=CommandLine.arguments
    try check((14...16).contains(a.count),"usage: hybrid_pilot LIB SHADER OUT KNOWN BACKEND D FAMILY BATCH BATCHES WORKERS SEED CAP MILLIS [MODE [COMPACT]]")
    let d=Int(a[6])!,family=Int(a[7])!,n=Int(a[8])!,batches=Int(a[9])!,workers=Int(a[10])!,seed=UInt64(a[11])!,cap=UInt32(a[12])!,millis=UInt64(a[13])!
    try check(["cpu","gpu"].contains(a[5]) && (2...10).contains(workers) && (1...1048576).contains(n) && (1...128).contains(batches) && ((9...17).contains(d) && (0...1).contains(family) || d==18 && (2...3).contains(family) || d==25 && family==4),"invalid backend/budget/batch/family")
    let producers=workers>=8 ? 2 : 1,validators=workers-producers
    let mode=a.count>=15 ? Int(a[14])! : 1,compact=a.count==16 && a[15]=="1"
    try check(mode==1 || mode==2,"mode must be 1 (screen) or 2 (screen plus early braid witness)")
    guard let lib=dlopen(a[1],RTLD_NOW|RTLD_LOCAL) else {throw NSError(domain:"dlopen",code:1)}
    defer {dlclose(lib)}
    func symbol<T>(_ name:String,_ type:T.Type)throws->T {
        guard let p=dlsym(lib,name) else {throw NSError(domain:name,code:1)}
        return unsafeBitCast(p,to:type)
    }
    let open=try symbol("hybrid_open",OpenSession.self),close=try symbol("hybrid_close",CloseSession.self)
    let fill=try symbol("hybrid_fill",Fill.self),screen=try symbol(mode==2 ? "hybrid_cpu_probe" : "hybrid_cpu_screen",ScreenCPU.self),validate=try symbol("hybrid_validate",Validate.self)
    guard let session=a[3].withCString({out in a[4].withCString {known in open(out,known)}}) else {throw NSError(domain:"session",code:1)}
    defer {close(session)}
    guard let device=MTLCreateSystemDefaultDevice(),let metalQueue=device.makeCommandQueue() else {throw NSError(domain:"Metal",code:1)}
    let shader=try device.makeLibrary(source:String(contentsOfFile:a[2],encoding:.utf8),options:nil)
    let pipeline=try device.makeComputePipelineState(function:shader.makeFunction(name:"screen")!)
    if mode==2 {
        // Independent Python oracle, exact S4 states after all 64 braid moves.
        // Checked before any timed workload; not merely agreement on rejections.
        let root=URL(fileURLWithPath:a[1]).deletingLastPathComponent()
        let oracle=try Data(contentsOf:root.appendingPathComponent("ORACLE_WALK_INPUT.bin")),expected=try Data(contentsOf:root.appendingPathComponent("ORACLE_WALK_OUTPUT.bin"))
        let count=oracle.count/12
        try check(oracle.count==count*12 && expected.count==count*4,"invalid independent walk oracle sizes")
        let diagnostic=try symbol("bench_screen",Diagnostic.self)
        let input=device.makeBuffer(length:oracle.count,options:.storageModeShared)!,gpu=device.makeBuffer(length:expected.count,options:.storageModeShared)!,cpu=device.makeBuffer(length:expected.count,options:.storageModeShared)!
        _ = oracle.withUnsafeBytes {memcpy(input.contents(),$0.baseAddress!,$0.count)}
        try check(diagnostic(4,UInt32(count),3,input.contents(),cpu.contents())==0,"CPU walk oracle failed")
        let command=metalQueue.makeCommandBuffer()!,encoder=command.makeComputeCommandEncoder()!
        encoder.setComputePipelineState(pipeline);encoder.setBuffer(input,offset:0,index:0);encoder.setBuffer(gpu,offset:0,index:1)
        let params=[UInt32(4),UInt32(count),UInt32(3)]
        params.withUnsafeBytes {encoder.setBytes($0.baseAddress!,length:$0.count,index:2)}
        encoder.dispatchThreads(MTLSize(width:count,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:128,height:1,depth:1));encoder.endEncoding();command.commit();command.waitUntilCompleted()
        try check(command.status == .completed,"GPU walk oracle command failed")
        try expected.withUnsafeBytes {try check(memcmp($0.baseAddress!,cpu.contents(),$0.count)==0 && memcmp($0.baseAddress!,gpu.contents(),$0.count)==0,"independent exhaustive walk oracle mismatch")}
        print("independent 64-move S4 oracle passed: \(count) tuples")
    }
    let slots=3,free=DispatchSemaphore(value:slots),finished=DispatchGroup()
    let validationQueue=DispatchQueue(label:"bounded-validation")
    let lock=NSLock();var rows=[[String:Any]](),errors=[String](),available=Array(0..<slots)
    defer {finished.wait()}
    var inputs=[MTLBuffer](),outputs=[MTLBuffer]()
    for _ in 0..<slots {
        guard let input=device.makeBuffer(length:n*3*d,options:.storageModeShared),let output=device.makeBuffer(length:n*4,options:.storageModeShared) else {throw NSError(domain:"buffers",code:1)}
        inputs.append(input);outputs.append(output)
    }
    // Explicit free-slot ownership: callback order need not equal submission
    // order. A slot stays owned until CPU validation finishes.
    let start=now()
    for batch in 0..<batches {
        free.wait();lock.lock();let failed = !errors.isEmpty;lock.unlock()
        if failed {free.signal();break;}
        lock.lock();let slot=available.removeLast();lock.unlock()
        let input=inputs[slot],output=outputs[slot]
        let batchSeed=seed &+ UInt64(batch)*3571
        let begin=now();try check(fill(session,UInt32(d),UInt32(n),UInt32(family),batchSeed,UInt32(producers),input.contents())==0,"generation failed")
        let generated=now()-begin
        // Save before submitting expensive work. On interruption the complete
        // batch and deterministic seed are retained, even without a summary.
        let prefix=String(format:"%@/batch_%04d",a[3],batch)
        if compact {
            let bytes=Data(bytesNoCopy:input.contents(),count:n*3*d,deallocator:.none)
            let metadata:[String:Any] = ["format":"HURWITZ_GENERATED_INPUT_V1","degree":d,"family":family,"count":n,"seed":String(batchSeed),"sha256":digest(bytes)]
            try JSONSerialization.data(withJSONObject:metadata,options:[.sortedKeys]).write(to:URL(fileURLWithPath:prefix+".input.json"),options:.atomic)
        } else {try Data(bytes:input.contents(),count:n*3*d).write(to:URL(fileURLWithPath:prefix+".input.bin"),options:.atomic)}
        finished.enter()
        func consume(_ screened:Double,_ deviceTime:Double) {
            validationQueue.async {
                defer {lock.lock();available.append(slot);lock.unlock();free.signal();finished.leave()}
                do {
                    if compact {try packedFlags(output.contents().assumingMemoryBound(to:UInt32.self),n).write(to:URL(fileURLWithPath:prefix+".flags.bits"),options:.atomic)}
                    else {try Data(bytes:output.contents(),count:n*4).write(to:URL(fileURLWithPath:prefix+".flags.bin"),options:.atomic)}
                    var stats=[UInt64](repeating:0,count:12);let vstart=now()
                    let rc=stats.withUnsafeMutableBytes {validate(session,UInt32(d),UInt32(n),UInt32(validators),cap,millis,input.contents(),output.contents(),$0.baseAddress!)}
                    try check(rc==0,"validation failed")
                    let row:[String:Any] = ["batch":batch,"generated_seconds":generated,"screen_wall_seconds":screened,"gpu_device_seconds":deviceTime,"validation_seconds":now()-vstart,"stats":stats]
                    try JSONSerialization.data(withJSONObject:row,options:[.sortedKeys]).write(to:URL(fileURLWithPath:prefix+".json"),options:.atomic)
                    lock.lock();rows.append(row);lock.unlock()
                    print("batch=\(batch) survivors=\(stats[0]) new_generating=\(stats[5]) certified=\(stats[8]) unresolved=\(stats[9])")
                    fflush(stdout)
                } catch {lock.lock();errors.append(String(describing:error));lock.unlock()}
            }
        }
        if a[5]=="cpu" {
            // Screening and validation share the same worker budget/serial stage.
            validationQueue.async {
                let s=now();let rc=screen(UInt32(d),UInt32(n),UInt32(validators),input.contents(),output.contents())
                if rc != 0 {lock.lock();errors.append("CPU screen failed");available.append(slot);lock.unlock();free.signal();finished.leave();return}
                // Already on the validation queue: enqueuing consume maintains
                // submission order, but new CPU screens may sit ahead of it.
                consume(now()-s,0)
            }
        } else {
            guard let command=metalQueue.makeCommandBuffer(),let encoder=command.makeComputeCommandEncoder() else {lock.lock();errors.append("encoder creation failed");available.append(slot);lock.unlock();free.signal();finished.leave();break}
            encoder.setComputePipelineState(pipeline);encoder.setBuffer(input,offset:0,index:0);encoder.setBuffer(output,offset:0,index:1)
            let params=[UInt32(d),UInt32(n),UInt32(mode)]
            params.withUnsafeBytes {encoder.setBytes($0.baseAddress!,length:$0.count,index:2)}
            encoder.dispatchThreads(MTLSize(width:n,height:1,depth:1),threadsPerThreadgroup:MTLSize(width:min(128,pipeline.maxTotalThreadsPerThreadgroup),height:1,depth:1))
            encoder.endEncoding();let submitted=now()
            command.addCompletedHandler {cmd in
                if cmd.status != .completed {lock.lock();errors.append("GPU command failed: \(String(describing:cmd.error))");available.append(slot);lock.unlock();free.signal();finished.leave();return}
                consume(now()-submitted,cmd.gpuEndTime-cmd.gpuStartTime)
            }
            command.commit()
        }
    }
    finished.wait();let elapsed=now()-start
    try check(errors.isEmpty,errors.joined(separator:"; "))
    // Entire-batch exact CPU shadow audit AFTER timed pipeline; not a hidden
    // CPU workload inside the reported GPU timings. Saved flags audited too.
    let auditStart=now();var shadow=[UInt32](repeating:0,count:n)
    for batch in 0..<batches {
        let prefix=String(format:"%@/batch_%04d",a[3],batch)
        let bytes:Data,flags:Data
        if compact {
            let metadata=try JSONSerialization.jsonObject(with:Data(contentsOf:URL(fileURLWithPath:prefix+".input.json"))) as! [String:Any]
            try check(metadata["format"] as? String=="HURWITZ_GENERATED_INPUT_V1" && metadata["degree"] as? Int==d && metadata["family"] as? Int==family && metadata["count"] as? Int==n,"input recipe metadata mismatch")
            guard let savedSeed=UInt64(metadata["seed"] as? String ?? "") else {throw NSError(domain:"input-seed",code:1)}
            try check(savedSeed==seed &+ UInt64(batch)*3571,"input seed mismatch")
            try check(fill(session,UInt32(d),UInt32(n),UInt32(family),savedSeed,UInt32(workers),inputs[0].contents())==0,"input regeneration failed")
            bytes=Data(bytesNoCopy:inputs[0].contents(),count:n*3*d,deallocator:.none)
            try check(digest(bytes)==metadata["sha256"] as? String,"regenerated input SHA256 mismatch")
            flags=try Data(contentsOf:URL(fileURLWithPath:prefix+".flags.bits"))
        } else {
            bytes=try Data(contentsOf:URL(fileURLWithPath:prefix+".input.bin"))
            flags=try Data(contentsOf:URL(fileURLWithPath:prefix+".flags.bin"))
        }
        let rc=bytes.withUnsafeBytes {input in shadow.withUnsafeMutableBytes {out in screen(UInt32(d),UInt32(n),UInt32(workers),input.baseAddress!,out.baseAddress!)}}
        try check(rc==0,"shadow screen failed")
        let equal=compact ? flags.count==(n+7)/8 && (0..<n).allSatisfy({shadow[$0]==UInt32((flags[$0/8]>>UInt8($0%8))&1)}) : shadow.withUnsafeBytes {expected in flags.withUnsafeBytes {actual in actual.count==expected.count && memcmp(actual.baseAddress!,expected.baseAddress!,actual.count)==0}}
        try check(equal,"shadow mismatch")
    }
    let result:[String:Any] = ["backend":a[5],"screen_mode":mode,"input_storage":compact ? "frozen_generator_seed_sha256" : "raw","flags_storage":compact ? "packed_bits" : "raw","degree":d,"family":family,"batch_size":n,"batches":batches,"cpu_budget":workers,"producers":producers,"validators":validators,"slots":slots,"seed":seed,"pipeline_seconds":elapsed,"shadow_audit_seconds":now()-auditStart,"shadow_mismatches":0,"device":device.name,"rows":rows.sorted {($0["batch"] as! Int)<($1["batch"] as! Int)}]
    try JSONSerialization.data(withJSONObject:result,options:[.prettyPrinted,.sortedKeys]).write(to:URL(fileURLWithPath:a[3]+"/SUMMARY.json"),options:.atomic)
    print("COMPLETE backend=\(a[5]) cpu_budget=\(workers) seconds=\(elapsed)")
}
do {try run()} catch {fputs("\(error)\n",stderr);exit(1)}
