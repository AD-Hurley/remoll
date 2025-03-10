#include <iostream>
#include <vector>
#include <algorithm>

struct DataGrid{
  vector<double> Q2, W, cosThetaK;
};

int get_lower_index(vector<double> &x, double xP){
  if(xP <= x[0])
    return 0;
  else if(xP >= x.back())
    return x.size()-2;
  else
    return lower_bound(x.begin(), x.end(), xP) - x.begin() - 1;     // lower_bound gives values >= xP  
}

double InterpolateRespFn(double Q2, double W, double cosThetaK){
  return 0;
}

void linearInterpTest(){
  
  DataGrid LambdaDataGrid; //(x,y,z) coordinates for the real data points
  LambdaDataGrid.Q2.push_back(0);
  //LambdaDataGrid.Q2.push_back(0.5);
  LambdaDataGrid.Q2.push_back(1.8);
  LambdaDataGrid.W.push_back(1.63);
  //LambdaDataGrid.W.push_back(1.90);
  LambdaDataGrid.W.push_back(2.569);
  LambdaDataGrid.cosThetaK.push_back(-0.8);
  //LambdaDataGrid.cosThetaK.push_back(0.0);
  LambdaDataGrid.cosThetaK.push_back(0.9);

  double testData[2][2][2]; //value of cross section at each (Q2,W,cosThetaK) point
  testData[0][0][0] = 0;
  testData[0][0][1] = 10;
  testData[0][1][0] = 0;
  testData[0][1][1] = 1;
  testData[1][0][0] = 0;
  testData[1][0][1] = 2;
  testData[1][1][0] = 0;
  testData[1][1][1] = .1;
  
  //Test Thrown event (Q2,W,cosThetaK)
  double TestThrown[3]{1.78,1.65,0.85};

  int iQ2 = get_lower_index(LambdaDataGrid.Q2, TestThrown[0]);
  int iW = get_lower_index(LambdaDataGrid.W, TestThrown[1]);
  int icosK = get_lower_index(LambdaDataGrid.cosThetaK, TestThrown[2]);

	int jQ2 = iQ2+1;
	int jW = iW+1;
	int jcosK = icosK+1;

  cout << iQ2 << "," << iW << "," << icosK << endl;
	cout << jQ2 << "," << jW << "," << jcosK << endl;
  
  double wQ2a = (TestThrown[0] - LambdaDataGrid.Q2[iQ2])/(LambdaDataGrid.Q2[jQ2]-LambdaDataGrid.Q2[iQ2]);
  double wWa = (TestThrown[1] - LambdaDataGrid.W[iW])/(LambdaDataGrid.W[jW]-LambdaDataGrid.W[iW]);
  double wcosKa = (TestThrown[2] - LambdaDataGrid.cosThetaK[icosK])/(LambdaDataGrid.cosThetaK[jcosK]-LambdaDataGrid.cosThetaK[icosK]); 
  
  cout << endl;
  cout << wQ2a << "," << wWa << "," << wcosKa << endl;
  
  double wQ2b = 1 - wQ2a;
  double wWb = 1 - wWa;
  double wcosKb = 1 - wcosKa; 
  
  cout << endl;
  cout << wQ2b << "," << wWb << "," << wcosKb << endl;
 
 	double crossSection = ((testData[iQ2][iW][icosK]*wQ2b +                      // return testData[iQ2][iW][icosK]*wQ2b*wWb*wcosKb +
             						  testData[jQ2][iW][icosK]*wQ2a) * wWb +               //        testData[jQ2][iW][icosK]*wQ2a*wWb*wcosKb +
            						 (testData[iQ2][jW][icosK]*wQ2b +                      //        testData[iQ2][jW][icosK]*wQ2b*wWa*wcosKb +
             							testData[jQ2][jW][icosK]*wQ2a) * wWa) * wcosKb +     //        testData[jQ2][jW][icosK]*wQ2a*wWa*wcosKb +
           							((testData[iQ2][iW][jcosK]*wQ2b +                      //        testData[iQ2][iW][jcosK]*wQ2b*wWb*wcosKa +
             							testData[jQ2][iW][jcosK]*wQ2a) * wWb +               //        testData[jQ2][iW][jcosK]*wQ2a*wWb*wcosKa +    
            						 (testData[iQ2][jW][jcosK]*wQ2b +                      //        testData[iQ2][jW][jcosK]*wQ2b*wWa*wcosKa +
             							testData[jQ2][jW][jcosK]*wQ2a) * wWa) * wcosKa;      //        testData[jQ2][jW][jcosK]*wQ2a*wWa*wcosKa;
             
  cout << crossSection << endl;
}

