"""Recover the complete cone certificate printed in Appendix A."""
import base64
import hashlib
import json
from pathlib import Path
import zlib

DATA = """
eNrtmH1oVWUcx/k+5+055968k+dZ98bddL7AOcd7h7vOVNrdG6FtErbNCYa1CRK7YtQfE0pMnGQs
yXKl4cDEhdJGDCZmkIhtQTKDSHvDBRMXKQYVYUwyJtjvOXei87U/xEZdDuc8b+d5ns/9Ps/vd3/P
+dzEKNMQPRc+GE7Y2BPusOFKlF2GrkODnueGnu7iSK4J2bAQ7VXVsaTrQUgqRzo5EgWIQsB5h0Uh
40hggY03ONUkEFL1r7E5DmzMDKGhEKEhGgRzEXNsVGIV2hDNjyMqVSeBUvxuYlYwf4pXepAyyTWi
kB6Ns7hAAcCWSPm/bTXlJfsvVoQP7VqIHoYm+21bCg/7GH/EscvFsAnGDAO6oTGmYwRLaEJDA9dp
QGMEazosJM4r5DF0WjT3QsW6DSG/JsCdEUY3JHZbWB2gitVrY7I04JT7LB1CA2go5lM3AXdRDC5e


1jmNlAZ1KCuYI+F4ckvUcZxfxkmbuhmWw3E6JpCeNUl/MMLVleTnAlKdk9Bc+xFlavwsa+JryOhN
rPPhzCjCzixpc6ipeW0zMRZ2WaQ91xmg6WD1KJudpC4vhXUQIQngYRotc0mRGsdx6Tm2xTygfWdH
e3BE6+U4ytDOpxATpphpjGmJfstgYcyDm6ZljaxaplZFVOGVyMe8YRNhZYo3QVT08V0k3fQMKsN4
Eb185jFkamlDFOKwZYbfxGy7FemwRmy6RQipkqkRAvGw65JxhV2w1x1DN3uXY7+ODr5D+54zU+5E
foTYyx9HXjXEPH87hk227AU/5cL1NYRF9WFs50/AeRQ1iWsUzy70a4owbUZzIUEcxjAxL0fNr2aW
Yf2tCP7VUbNOyjTSAmkpxqy3HpIfoXIvsINXQHzKRjhp4dXBi1R5dXUy7buy/oJJ6xbRk3C96T6H
t8T2I4hAK7YdZ54st6Eu2gNHackqVGEBWjPF7ZBfkXU0JbGR7C2RIQNRxlIws1bWRjMxWRiwnjFh
KPlhmSo5BDuCLoLWgpvIyVxKBeJETznHw1U7x/9v8rd/gm0Vg1qPg0iMnGO6hMfIRfriS/sMljr9
dg+5OjSswPQfrMbWbnbEAW3CqPdcnFd5AxwbsBjhlcIrIR9nTRWhlHfRaly/GR+gnnzeWU7OKLEC
DQ9XYi+50jY8s9QhNxrDNhv1RhW5YzkLVp4MpZS7GPrTmAwU9iSAaLug4Wp46FWWW5ZrFKNmbnNM
IoiDr2Ok/Mhj5LlSridlnhoTBRQfuP5t6o7z1r2slYKojc/TLOp52sIgToAen5l56s9d1xrJIwZd
Emkp4d22dtjUKNyofgpxcpncknUS6VtqdqiAi2K6KFrojsvxzMSa3VYQ67XQM4SERDsmlqk9QWPG
BPnVbJfMzRWnzABsJQV+gqaugpsvFt1adWLUqCc9ylX05HrXs8e5vMJkIATEkzRzNj2r6FtqED2I
8cyAycHJgevFFD0Eq0dXHNIlQdw7Nw2SxVA4phlaElLpQ4GvLerSt609EGhWij1AsuV6NlBpJ3Et
eg/Xcr1WEDnPF2gkHZKBpDcU+4mWxqTbo5+qJuABleu7ivZOTYNmTqR7itR3SMPJ+fTXgBNtbDad
Z2QdhE/Bxw35/VbsMov1UDRehBayMJUMjYMvD04rQu3lieUB2sqGxnWDu8qELRVw00NKF5LOCXdr
HDQNFVZTiMGT+bRkMp5VWIUnwRvePdp7SL7zxLM5m/yk1MxU4ltk0z5VLlDqrYM6p22eWBwwg6MB
uHIWlHqemxVNqSaCM5iA9w/eITPNKXp/Fb1o5hS9v4rmBM0JOqkF7d+K90tPM3Tb6GNDc7/Qfrbc
DeNXdzL4SCfUp6K7pt+YTlouuenuU9+CqFmW3SPtvGQ88Ckf+I9UkcCp6s4t7P+g7h/mf39BH/yM
ZKpdpThZjb8BAF4yjg==
"""
EXPECTED = (
    "6c1b022e9351a42cc36343600484a7fd7bec93eb42065eb1f454580db5826d60"
)
stream = zlib.decompress(base64.b64decode(DATA))
records = []
previous = value = shift = 0
for byte in stream:
    value |= (byte & 127) << shift
    if byte & 128:
        shift += 7
        continue
    previous += value
    n = previous
    n, g = divmod(n, 2)
    n, z = divmod(n, 20)
    n, y = divmod(n, 20)
    n, x = divmod(n, 20)
    n, a = divmod(n, 7)
    tag, b = divmod(n, 5)
    records.append(["JSTV"[tag], [b, [0, 1, 2, 3, 4, 5, 7][a]],
                     [x, y, z], g])
    value = shift = 0
raw = (json.dumps(records) + "\n").encode("ascii")
if shift or len(records) != 5485 or hashlib.sha256(raw).hexdigest() != EXPECTED:
    raise ValueError("Certificate decoding or checksum failed.")
Path(__file__).with_name("original_cone_stratified.json").write_bytes(raw)
print("Recovered", len(records), "entries; SHA256:", EXPECTED)
