std::vector<std::string> fMap;
std::vector<std::string> fXAxis;
std::vector<std::string> fYAxis;
std::vector<int>         fNBins;
std::vector<double>      fXMin;
std::vector<double>      fXMax;

double *fAxis = new double[24];
fAxis =  {0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,250,300,350};
void loadfMap() {
  fMap.push_back("m_sv");     fNBins.push_back(23);   fXMin.push_back(0);    fXMax.push_back(300); fXAxis.push_back("m_{#tau#tau}  (GeV/c^{2})");  fYAxis.push_back("Events/10 GeV/c^{2}"); //fXArr.push_back(fAxis);
  fMap.push_back("m_vis");    fNBins.push_back(23);   fXMin.push_back(0);    fXMax.push_back(200); fXAxis.push_back("m_{vis}  (GeV/c^{2})");       fYAxis.push_back("Events/7.5 GeV/c^{2}"); //fXArr.push_back(fAxis);
}

int getId(std::string iStr) { 
  for(unsigned int i0 = 0; i0 < fMap.size(); i0++) { 
    if(iStr == fMap[i0]) return i0;
  }
  return 0;
}
const char*       getXAxis(std::string iStr) { return fXAxis[getId(iStr)].c_str();}
const char*       getYAxis(std::string iStr) { return fYAxis[getId(iStr)].c_str();}
int         getNBins(std::string iStr) { return fNBins[getId(iStr)];}
double      getXMin (std::string iStr) { return fXMin [getId(iStr)];}
double      getXMax (std::string iStr) { return fXMax [getId(iStr)];}
