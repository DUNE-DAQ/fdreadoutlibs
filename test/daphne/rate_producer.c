#define _POSIX_C_SOURCE 200809L
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <time.h>
#include <unistd.h>
#include <signal.h>
static volatile sig_atomic_t stopping;
static void stop(int signal_number){(void)signal_number;stopping=1;}
static uint64_t now(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return (uint64_t)t.tv_sec*1000000000+t.tv_nsec;}
int main(int argc,char **argv){
 if(argc!=3)return 2;
 unsigned long hz=strtoul(argv[1],0,10),n=strtoul(argv[2],0,10);
 if(n<2||n>2000000000||hz>1000000)return 2;
 int fd=open("/dev/mem",O_RDWR|O_SYNC);if(fd<0)return 3;
 volatile uint32_t *f=mmap(0,4096,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0x88000000);
 volatile uint32_t *s=mmap(0,4096,PROT_READ,MAP_SHARED,fd,0x94000000);
 volatile uint32_t *c=mmap(0,4096,PROT_READ,MAP_SHARED,fd,0xa0010000);
 if(f==MAP_FAILED||s==MAP_FAILED||c==MAP_FAILED)return 4;
 if(s[60]!=0x44415048||s[61]!=0x20000||s[62]!=1||s[63]!=0x14f56c3||f[13]!=0)return 5;
 for(int ch=0;ch<32;ch++)if(c[ch*8+7]&0x80000000u)return 6;
 signal(SIGTERM,stop);signal(SIGINT,stop);signal(SIGALRM,stop);alarm(120);
 unsigned long sent=0;
 uint64_t start=now(),previous=0,min_interval=UINT64_MAX,max_interval=0;
 if(hz){for(unsigned long i=0;i<n&&!stopping;i++){
  uint64_t deadline=i?previous+1000000000ull/hz:start;
  uint64_t stamp;while((stamp=now())<deadline&&!stopping){}
  if(stopping)break;
  f[2]=0xbaba;sent++;
  if(i){uint64_t interval=stamp-previous;if(interval<min_interval)min_interval=interval;if(interval>max_interval)max_interval=interval;}
  previous=stamp;
 }}
 else {for(unsigned long i=0;i<n&&!stopping;i++){f[2]=0xbaba;sent++;}}
 __asm__ volatile("dsb sy" ::: "memory");
 uint64_t end=now();
 printf("{\"target_hz\":%lu,\"writes\":%lu,\"elapsed_ns\":%llu,\"write_rate_hz\":%.3f,\"min_interval_ns\":%llu,\"max_interval_ns\":%llu}\n",hz,sent,(unsigned long long)(end-start),(double)(hz?(sent?sent-1:0):sent)*1e9/(end-start),(unsigned long long)(hz?min_interval:0),(unsigned long long)max_interval);
 munmap((void*)f,4096);munmap((void*)s,4096);munmap((void*)c,4096);close(fd);return 0;
}
