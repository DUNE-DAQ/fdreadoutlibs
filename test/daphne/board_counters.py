import mmap,os,struct,subprocess,json,time
fd=os.open('/dev/mem',os.O_RDONLY|os.O_SYNC)
bank=mmap.mmap(fd,4096,access=mmap.ACCESS_READ,offset=0xa0010000)
def count(o):
 for _ in range(10):
  hi=struct.unpack_from('<I',bank,o+4)[0];lo=struct.unpack_from('<I',bank,o)[0]
  if hi==struct.unpack_from('<I',bank,o+4)[0]: return hi<<32|lo
 raise RuntimeError('incoherent counter')
def snapshot():
 return [dict(channel=ch,records=count(ch*32+4),busy=count(ch*32+12),full=count(ch*32+20),continuations=count(2048+ch*32),continuation_rejected=count(2056+ch*32),merged=count(2064+ch*32),overflow=count(2072+ch*32)) for ch in range(32)]

print(json.dumps({"unix":time.time(),"monotonic":time.monotonic(),"channels":snapshot()}))
bank.close();os.close(fd)
