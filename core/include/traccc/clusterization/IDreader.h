#ifndef IDreader_h
#define IDreader_h

#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iostream>

using namespace std;

class IDreader {
public:

  IDreader(string filename);
  ~IDreader();
  unsigned long int GiveDetectorID(int lineID);
  unsigned long int Give50muPixelDetectorID(int lineID, int ifpixel=1);

  int NIDs(){return m_Nid;};

private:

  int m_Nid;
  vector<unsigned long int> m_GeoID;
  vector<unsigned long int> m_50muGeoID;  //ID of 50micron x 50micron pixels
  vector<unsigned long int> m_elseGeoID;

  const int min50muVolumeID=7;
  const int max50muVolumeID=9;

  const int m_NactiveVolumeID=9;
  const int m_activeVolumeID[9]={7,8,9,12,13,14,16,17,18};
  //traccc volume IDs for which config is given in config file. 
};
#endif

IDreader::IDreader(string filename){

  cout<<"IDreader will read traccc detector csv file: "<<filename<<endl;
  
  ifstream ifs;
  ifs.open(filename.data());

  string line;
  string separator=",";
  int iline=0;
  m_GeoID.clear();
  m_50muGeoID.clear();
  m_elseGeoID.clear();
  while(getline(ifs,line)){
    iline++;
    if(iline==1) continue;   //skip header line

    int pos1 = line.find(separator);
    string geoid = line.substr(0, pos1);
    m_GeoID.push_back(stoul(geoid));
    int pos2 = line.find(separator,pos1+1);
    string volumeid = line.substr(pos1+1,pos2-pos1-1);

    bool ifactive=false;
    for(int i=0; i<m_NactiveVolumeID; i++){
      if(stoi(volumeid)==m_activeVolumeID[i]){
	ifactive=true;
	break;
      }
    }
    if(!ifactive) continue;
    
    if(stoi(volumeid)>=min50muVolumeID && stoi(volumeid)<=max50muVolumeID)
      m_50muGeoID.push_back(stoul(geoid));
    else
      m_elseGeoID.push_back(stoul(geoid));
    
  }

  m_Nid = m_GeoID.size();
  cout<<"IDreader: "<<m_Nid<<" IDs are read."<<endl;
  cout<<"IDreader: "<<(m_50muGeoID.size()+m_elseGeoID.size())
      <<" IDs are considered to be active."<<endl;
}

IDreader::~IDreader(){

  m_GeoID.clear();
}


unsigned long int IDreader::GiveDetectorID(int lineID){

  if(lineID >= m_Nid){
    
    cout<<"line ID is too large."<<endl;
    return -1;

  }

  return m_GeoID[lineID];

}
  
unsigned long int IDreader::Give50muPixelDetectorID(int lineID, int ifpixel){

  int key=lineID;

  if(ifpixel==1){
    if(key >= m_50muGeoID.size()){
      key=lineID%m_50muGeoID.size();
      
    }
    return m_50muGeoID[key];

    
  }else{
    if(key >= m_elseGeoID.size()){
      key=lineID%m_elseGeoID.size();
    }
    return m_elseGeoID[key];

  }


}


