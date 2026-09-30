#include "fdreadoutlibs/pds/ActivityBuilder.hpp"
#include <fstream>
#include <iostream>
#include <chrono>
#include <vector>
#include <tuple>
using namespace dunedaq::fdreadoutlibs::pds;
int main(int argc,char**argv){
 if(argc!=3)return 2;
 try{
  std::ifstream input(argv[1],std::ios::binary|std::ios::ate);if(!input)throw std::runtime_error("Cannot open input");
  auto bytes=input.tellg();if(bytes<=0||uint64_t(bytes)%sizeof(TP))throw std::runtime_error("Bad TP binary size");
  std::vector<TP> tps(uint64_t(bytes)/sizeof(TP));input.seekg(0);input.read(reinterpret_cast<char*>(tps.data()),bytes);if(!input)throw std::runtime_error("Read failed");
  auto a=std::chrono::steady_clock::now();
  std::sort(tps.begin(),tps.end(),[](const TP&x,const TP&y){return std::tie(x.time_start,x.channel)<std::tie(y.time_start,y.channel);});
  for(size_t i=1;i<tps.size();++i)if(tps[i].time_start==tps[i-1].time_start&&tps[i].channel==tps[i-1].channel)throw std::runtime_error("Duplicate TP channel/start");
  auto b=std::chrono::steady_clock::now();
  CoincidenceBuilder coincidence({625,6,100000,1});PromptLightBuilder prompt({7,625,1});
  std::ofstream output(argv[2]);uint64_t activities=0,light_windows=0,total_integral=0;size_t max_inputs=0;
  for(const auto&tp:tps){
   if(auto ta=coincidence.push(tp)){
    ++activities;max_inputs=std::max(max_inputs,ta->inputs.size());
    output<<"{\"kind\":\"ta\",\"time_start\":"<<ta->time_start<<",\"time_end\":"<<ta->time_end<<",\"inputs\":"<<ta->inputs.size()<<",\"integral\":"<<ta->adc_integral<<"}\n";
   }
   if(auto light=prompt.push(tp)){++light_windows;total_integral+=light->total_integral;}
  }
  if(auto light=prompt.finish()){++light_windows;total_integral+=light->total_integral;}
  auto c=std::chrono::steady_clock::now();double sort_s=std::chrono::duration<double>(b-a).count(),build_s=std::chrono::duration<double>(c-b).count();
  std::cout<<"{\"input_tps\":"<<tps.size()<<",\"tp_size\":"<<sizeof(TP)<<",\"activities\":"<<activities<<",\"light_windows\":"<<light_windows<<",\"total_integral\":"<<total_integral<<",\"max_ta_inputs\":"<<max_inputs<<",\"sort_seconds\":"<<sort_s<<",\"builder_seconds\":"<<build_s<<",\"builder_tps_per_second\":"<<tps.size()/build_s<<",\"sort_and_build_tps_per_second\":"<<tps.size()/(sort_s+build_s)<<",\"minimum_channels\":6,\"window_us\":10,\"offline_sorted_replay\":true}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
}
