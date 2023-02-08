
// http://egloos.zum.com/nicepan/v/1157133

#include <iostream>
#include <string.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <typeinfo>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <complex.h> 

// 0 : idx
// 1 : act_tor
// 2 : tar_tor
// 3 : act_pos
// 4 : tar_pos
// 5 : diff_pos
// 6 : tar_vel
// 7 : act_vel
// 8 : calc_vel
// 9 : tar_accel
// 10 : tar_power(W)
// 11 : act_poser(W)

using namespace std;

 double string_to_double( const string& s )
 {
   istringstream i(s);
   double x;
   if (!(i >> x))
     return 0;
   return x;
 } 

int main(int argc, char** argv) {

	double Fs = 4000.0;     // Sampling frequency 4[kHz]
	double Fn = Fs/2;       // Nyquist  frequency 2[kHz]
	double Et = 10.0;       // End time
	int L = (int) Et*Fs;    // Length of signal
 	const int nbin = 65536; // Next power of 2 from length of signal

    int idx_cnt, at_cnt, tt_cnt, ap_cnt, tp_cnt, tv_cnt, av_cnt = 0;

	double idx[nbin]={};
	double act_tor[nbin]={};
	double tar_tor[nbin]={};
	double act_pos[nbin]={};
	double tar_pos[nbin]={};
	double tar_vel[nbin]={};
	double act_vel[nbin]={};

	double data_act[nbin]={};
	double data_tar[nbin]={};

	int IndexNum = (int)(nbin/2);
	double FrqVec[IndexNum+1]={};

	for (int i = 0; i < IndexNum + 1; ++i ){
		FrqVec[i] = (i*Fn)/(IndexNum);
	}

	double TF_Mag[nbin], TF_Phase[nbin];
	double act_Mag, tar_Mag;
	double act_Phase, tar_Phase;

	static const double PI = 3.141592653589793;

    fstream fs;
    string str_buf;         
    string data;
    vector<string> DATA; 
    vector<string>::iterator vi; 

    fs.open("pos_200hz.csv",ios::in);

    for (int i = 0; i < L+1; ++i ) // Contain first information index (L+1)
    {
	    getline(fs,str_buf);
	    stringstream sstream(str_buf);

	    for (int j = 0; j < 12; ++j )
	    	{
	    		getline(sstream,data,',');
		   		DATA.push_back(data);
	    	}
    }
	
	for (int k = 0; k < DATA.size(); k++)
	{	
		if ( k > 11)
		{
			if (k%12 == 0){
				idx[idx_cnt] = string_to_double(DATA[k]);
				idx_cnt = idx_cnt + 1;
			}
			if (k%12 == 1){
				act_tor[at_cnt] = string_to_double(DATA[k]);
				at_cnt = at_cnt + 1;
			}
			if (k%12 == 2){
				tar_tor[tt_cnt] = string_to_double(DATA[k]);
				tt_cnt = tt_cnt + 1;
			}
			if (k%12 == 3){
				act_pos[ap_cnt] = string_to_double(DATA[k]);
				ap_cnt = ap_cnt + 1;
			}
			if (k%12 == 4){
				tar_pos[tp_cnt] = string_to_double(DATA[k]);
				tp_cnt = tp_cnt + 1;
			}
			if (k%12 == 6){
				tar_vel[tv_cnt] = string_to_double(DATA[k]);
				tv_cnt = tv_cnt + 1;
			}
			if (k%12 == 7){
				act_vel[av_cnt] = string_to_double(DATA[k]);
				av_cnt = av_cnt + 1;
			}
		}
	} 

    fs.close();

    fftw_complex *act_ini, *act_fin;
    act_ini = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbin);
    act_fin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbin);

    fftw_complex *tar_ini, *tar_fin;
    tar_ini = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbin);
    tar_fin = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nbin);

	/* 푸리에 변환을 위한 fftw_plan형 변수를 선언하고 초기화 */
	fftw_plan act_fft;
	act_fft = fftw_plan_dft_1d(nbin, act_ini, act_fin, FFTW_FORWARD, FFTW_ESTIMATE);

	fftw_plan tar_fft;
	tar_fft = fftw_plan_dft_1d(nbin, tar_ini, tar_fin, FFTW_FORWARD, FFTW_ESTIMATE);

	// Data select & copy
	memcpy(&data_act, &act_pos, sizeof(act_pos));
	memcpy(&data_tar, &tar_pos, sizeof(tar_pos));

	// 변환 대상이 되는 함수를 지정
	for (int iq = 0; iq < nbin; iq++) {
	    act_ini[iq][0] = data_act[iq];   // Real
	    act_ini[iq][1] = 0.0;            // Imag
	    tar_ini[iq][0] = data_tar[iq];   // Real
	    tar_ini[iq][1] = 0.0;            // Imag
	}

	// Do FFT
	fftw_execute(act_fft);
	fftw_execute(tar_fft);

	for (int n = 0; n < nbin ; n++){

		act_Mag   = sqrt(pow(act_fin[n][0], 2.0) + pow(act_fin[n][1], 2.0));
		act_Phase = atan2(act_fin[n][1], act_fin[n][0]);

		tar_Mag   = sqrt(pow(tar_fin[n][0], 2.0) + pow(tar_fin[n][1], 2.0));
		tar_Phase = atan2(tar_fin[n][1], tar_fin[n][0]);

		TF_Mag[n]   = 20*log10(act_Mag/tar_Mag);
		TF_Phase[n] = (act_Phase - tar_Phase)*180.0/PI;

		// cout << n << " TF_Mag : " << act_Mag/tar_Mag << " TF_Phase : " << act_Phase - tar_Phase << endl ;
	} 

	// 동적으로 할당된 메모리 해제
	fftw_destroy_plan(act_fft);
	fftw_destroy_plan(tar_fft);
	fftw_free(act_ini);
	fftw_free(act_fin);
	fftw_free(tar_ini);
	fftw_free(tar_fin);
    
    return 0;
}