"""Local QA only: reproduce the public read-only URL mappings."""
import argparse
from functools import partial
from http.server import SimpleHTTPRequestHandler,ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit


class Handler(SimpleHTTPRequestHandler):
    extensions_map={**SimpleHTTPRequestHandler.extensions_map,'.mjs':'application/javascript'}
    def do_GET(self):
        path=urlsplit(self.path).path
        if path in ['/all-smooth','/all-smooth/','/','/all-one-piece','/all-one-piece/']:
            self.path='/ui/index.html'
        elif path.startswith('/all-smooth-assets/'):
            self.path='/ui/'+path.rsplit('/',1)[1]
        return super().do_GET()


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--release',type=Path,required=True);p.add_argument('--port',type=int,default=8769)
    args=p.parse_args();server=ThreadingHTTPServer(('127.0.0.1',args.port),partial(Handler,directory=str(args.release.resolve())))
    print(f'Local catalog: http://127.0.0.1:{args.port}/all-smooth',flush=True);server.serve_forever()
