#!/usr/bin/env python3

import sys
import requests
import argparse

if __name__ != '__main__':
    sys.exit()

argparser = argparse.ArgumentParser()
argparser.add_argument('username')
argparser.add_argument('password')
argparser.add_argument('--user-agent')
args = argparser.parse_args()
if not args.user_agent:
    args.user_agent = 'its-client-cmd/1.0'

s = requests.Session()
s.headers['User-Agent'] = args.user_agent
s.get('https://its.pku.edu.cn') # set cookies
r = s.post('https://its.pku.edu.cn/cas/webLogin', data=dict(iprange='yes', username=args.username, password=args.password))
#import pdb; pdb.set_trace()
sys.exit(int(r.status_code != 200))


