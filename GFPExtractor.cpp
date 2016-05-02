/***************************************************************************
 
 *  Copyright 2015 Hadi Harb                                                *
 *                                                                          *
 *  Licensed under the Apache License, Version 2.0 (the "License");         *
 *  you may not use this file except in compliance with the License.        *
 *  You may obtain a copy of the License at                                 *
 *                                                                          *
 *  http://www.apache.org/licenses/LICENSE-2.0                              *
 *                                                                          *
 *  Unless required by applicable law or agreed to in writing, software     *
 *  distributed under the License is distributed on an "AS IS" BASIS,       *
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.*
 *  See the License for the specific language governing permissions and     *
 *  limitations under the License.                                          *
 
 ***************************************************************************/



#include <math.h>
#include <stdio.h>
#include <string.h>

#include "GFPExtractor.h"


// *********************************************************************************************************


CFPExtractor::CFPExtractor(int inWinSize, int inWinStep, int inIntWinSize, int inFracUsed, float inOverlap)
{
    m_win_size = inWinSize;
    m_win_step = inWinStep;
    m_intwin_size = inIntWinSize;
    m_frac_used = inFracUsed;
    m_intwin_overlap = (int)(inOverlap*m_intwin_size);
    m_cnt_xbuf = m_cnt_featbuf = m_cnt_resultbuf = 0;
    double pi=3.14159265359;
 
    
    if (m_win_size<=256)
    {
        m_FFT = new GFFT(256);
        m_win_size = 256;
    } else {
        m_FFT = new GFFT(512);
        m_win_size = 512;
    }
    for (int i=0; i<GF_WORKBUF_SIZE; i++)
	m_fftbuf[i] = new double[m_win_size];
    
// init ITW Hamming window
    m_ITWhamming = new float[m_intwin_size + 2*m_intwin_overlap];
    m_intwin_total = m_intwin_size + 2*m_intwin_overlap;
    for (int i=0; i<m_intwin_total; i++)
    {
	if (i>=m_intwin_overlap && i<m_intwin_size+m_intwin_overlap)
	    m_ITWhamming[i] = 1.0f;
	else
	    m_ITWhamming[i] = 0.54f - 0.46*cos(2*M_PI*(float)i / (float)(m_intwin_total-1));
    }
    
    for(int k=0;k<GF_PASSFREQ2;k++)
    {
        dfttC[k]=new float[m_intwin_total];
        dfttS[k]=new float[m_intwin_total];
        for(int n=0;n<m_intwin_total;n++)
        {
            dfttS[k][n]=sin(-n*k*2*pi/m_intwin_total);
            dfttC[k][n]=cos(-n*k*2*pi/m_intwin_total);
        }
    }
    
	
// init signal Hamming window
    m_hamming = new double[m_win_size];
    for (int i=0; i<m_win_size; i++)
	m_hamming[i] = 0.54 - 0.46*cos(2*M_PI*(double)i/(double)(m_win_size-1));
	
    m_tmp_buf = new double[m_win_size];
    m_sliding_norm = new double[GF_SLIDING_NORM_SIZE];
    m_sliding_norm1 = new double[GF_SLIDING_NORM_SIZE];
    for (int i=0; i<GF_SLIDING_NORM_SIZE; i++)
    {
	m_sliding_norm[i] = 0.0;
	m_sliding_norm1[i] = 0.0;
    }
    m_sliding_norm_ptr = 0;
    m_sliding_norm_over = false;

    for (int i=0; i<13; i++) m_filt[i]=0.003;
    for (int i=13; i<113; i++) m_filt[i]=1;
    for (int i=113; i<200; i++) m_filt[i]=0.003;
    int L=0;
    for (int i=0;i<GF_PASSFREQ/GF_PASSBND;i++)
    {
        for (int j=0;j<GF_PASSFREQ/GF_PASSBND;j++)
        {
            L=j%2;
            if(L>1)L=1;
            Gbr[0].mat[i][j]=L;
        }
    }
    for (int i=0;i<GF_PASSFREQ/GF_PASSBND;i++)
    {
        for (int j=0;j<GF_PASSFREQ/GF_PASSBND;j++)
        {
            L=j%4;
            if(L<2)L=0;
            else L=1;
            Gbr[1].mat[i][j]=L;
        }
    }
    for (int i=0;i<GF_PASSFREQ/GF_PASSBND;i++)
    {
        for (int j=0;j<GF_PASSFREQ/GF_PASSBND;j++)
        {
            L=j%6;
            if(L<3)L=1;
            else L=0;
            Gbr[2].mat[i][j]=L;
        }
    }
    for (int i=0;i<GF_PASSFREQ/GF_PASSBND;i++)
    {
        for (int j=0;j<GF_PASSFREQ/GF_PASSBND;j++)
        {
            L=i%2;
            if(L>1)L=1;
            Gbr[3].mat[i][j]=L;
        }
    }
    for (int i=0;i<GF_PASSFREQ/GF_PASSBND;i++)
    {
        for (int j=0;j<GF_PASSFREQ/GF_PASSBND;j++)
        {
            L=i%4;
            if(L<2)L=0;
            else L=1;
            Gbr[4].mat[i][j]=L;
        }
    }
    for (int i=0;i<GF_PASSFREQ/GF_PASSBND;i++)
    {
        for (int j=0;j<GF_PASSFREQ/GF_PASSBND;j++)
        {
            L=i%6;
            if(L<3)L=0;
            else L=1;
            Gbr[5].mat[i][j]=L;
        }
    }
    for (int i=0;i<GF_PASSFREQ/GF_PASSBND;i++)
    {
        for (int j=0;j<GF_PASSFREQ/GF_PASSBND;j++)
        {
            L=(i+j)%2;
            if(L>1)L=1;
            Gbr[6].mat[i][j]=L;
        }
    }
    for (int i=0;i<GF_PASSFREQ/GF_PASSBND;i++)
    {
        for (int j=0;j<GF_PASSFREQ/GF_PASSBND;j++)
        {
            L=(i+j)%4;
            if(L<2)L=0;
            else L=1;
            Gbr[7].mat[i][j]=L;
        }
    }
}

// *********************************************************************************************************

CFPExtractor::~CFPExtractor()
{
    delete m_sliding_norm1;
    delete m_sliding_norm;
    delete m_tmp_buf;
    delete m_hamming;
    delete m_ITWhamming;
    for (int i=0; i<GF_WORKBUF_SIZE; i++)
	delete m_fftbuf[i];
    delete m_FFT;
}

// *********************************************************************************************************

bool CFPExtractor::MakeFeatureVectorFromX(TXVector* inBuf, int inBufSize, int inPos, TFeatureVector* out)
{
    if (inBufSize-inPos < m_intwin_total || !out || !inBuf) return false;
    
    float mean, var;
    for (int k=0; k<GF_PASSFREQ; k++)
    {
	mean = var = 0.0f;
        for (int l=inPos; l<inPos+m_intwin_total; l+=1)
        {
            float tmp = inBuf[l][k]*m_ITWhamming[l-inPos];
            mean += tmp;
            var += tmp*tmp;
        }
        out->mean[k] = mean / m_intwin_size;
        out->var[k] = (var / m_intwin_size) - pow (out->mean[k],2);
    }

// pitch

    float maxm=0.0f,maxm1=0.0f,maxm2=0.0f;
    float maxK=0.0f,maxK1=0.0f,maxK2=0.0f;
    for (int k=0; k<GF_PASSFREQ; k++)
    {
        float tmp = out->mean[k]*(k+1);
        maxm += tmp;
        maxK += out->mean[k];
        if (k>=1 && k<GF_PASSFREQ/2)
        {
	    maxm1 += tmp;
	    maxK1 += out->mean[k];
        } else if (k>=GF_PASSFREQ/2)
        {
	    maxm2 += tmp;
	    maxK2 += out->mean[k];
        }
    }
//    maxK += 0.0001;
//    maxK1 += 0.0001;
//    maxK2 += 0.0001;
    out->maxmean.value = maxm / (maxK+0.01);
    out->maxmean.index = (int)maxK;
    out->maxmeanLow.value = maxm1 / (maxK1+0.01);
    out->maxmeanLow.index = (int)maxK1;
    out->maxmeanHi.value= maxm2 / (maxK2+0.01);
    out->maxmeanHi.index = (int)maxK2;
    
    return true;
}

// *********************************************************************************************************

bool CFPExtractor::MakeFeatureVector2FromX(TXVector* inBuf, int inBufSize, int inPos, TFeatureVector2* out)
{
    if (inBufSize-inPos < m_intwin_total || !out || !inBuf) return false;
    
    double pin;
    float buffers, bufferc, buffers1, bufferc1, bufferx, bufferz, bufferI;
    bufferc=buffers=bufferc1=buffers1=bufferx=bufferz=bufferI=0;
    pin=3.14159265359;
    
    for(int k=0; k<GF_PASSFREQ; k+=GF_PASSFREQ/GF_PASSBND)
    {
        for(int i=0; i<GF_PASSFREQ2; i++)
        {
            bufferc=buffers=0;
            for (int j=inPos; j<inPos+m_intwin_total; j+=1)
            {
                bufferx = 0;
                for(int l=k;l<k+GF_PASSFREQ/GF_PASSBND;l++)
                {
                    bufferx+=inBuf[j][l];
                }
                bufferx/=GF_PASSFREQ/GF_PASSBND;
                buffers=buffers+bufferx*dfttS[i][j-inPos];
                bufferc=bufferc+bufferx*dfttC[i][j-inPos];
            }
            out->feat[k/(GF_PASSFREQ/GF_PASSBND)][i]=sqrt(buffers*buffers+bufferc*bufferc);
        }
    }
 /*
    for(int k=0; k<GF_PASSFREQ; k+=GF_PASSFREQ/GF_PASSBND)
    {
        for (int j=inPos; j<inPos+m_intwin_total; j+=1)
        {
            bufferx = 0;
            bufferz = 0;
            bufferI = 0;
            for(int l=k;l<k+GF_PASSFREQ/GF_PASSBND;l++)
            {
                bufferx+=inBuf[j][l];
                bufferI+=(l-k+1)*inBuf[j][l];
            }
            bufferz=bufferI/(0.01+bufferx);
            buffers1+=bufferz;
        }
        out->feat2[k/(GF_PASSFREQ/GF_PASSBND)]=buffers1/m_intwin_total;
        buffers1 = 0;
    }
    */
    /*
    float tmp = 0, tmp2 = 0, tmp3 = 0;
    for(int k=0; k<GF_PASSFREQ; k+=GF_PASSFREQ/GF_PASSBND)
    {
        tmp3=0;
        for (int j=inPos; j<inPos+m_intwin_total; j+=1)
        {
            tmp = tmp2 = 0;
            for(int l=k;l<k+GF_PASSFREQ/GF_PASSBND;l++){ tmp+=inBuf[j][l];tmp2+=(l-k+1)*inBuf[j][l];}
            tmp3+=tmp2/(0.01+tmp);
        }
        tmp3/=m_intwin_total;
        out->feat4[k/(GF_PASSFREQ/GF_PASSBND)]=tmp3;
    }
   
    for(int k=0; k<GF_PASSFREQ; k+=2*GF_PASSFREQ/GF_PASSBND)
    {
        for(int i=0; i<GF_PASSFREQ2; i++)
        {
            bufferc=buffers=0;
            for (int j=inPos; j<inPos+m_intwin_total; j+=1)
            {
                bufferx = 0;
                for(int l=k;l<k+2*GF_PASSFREQ/GF_PASSBND;l++) bufferx+=inBuf[j][l];
                bufferx/=2*GF_PASSFREQ/GF_PASSBND;
                buffers=buffers+bufferx*dfttS[i][j-inPos];
                bufferc=bufferc+bufferx*dfttC[i][j-inPos];
                
            }
            out->feat2[k/(2*GF_PASSFREQ/GF_PASSBND)][i]=sqrt(buffers*buffers+bufferc*bufferc);
        }
    }
    
    for(int k=0; k<GF_PASSFREQ; k+=4*GF_PASSFREQ/GF_PASSBND)
    {
        for(int i=0; i<GF_PASSFREQ2; i++)
        {
            bufferc=buffers=0;
            for (int j=inPos; j<inPos+m_intwin_total; j+=1)
            {
                bufferx = 0;
                for(int l=k;l<k+4*GF_PASSFREQ/GF_PASSBND;l++) bufferx+=inBuf[j][l];
                bufferx/=4*GF_PASSFREQ/GF_PASSBND;
                buffers=buffers+bufferx*dfttS[i][j-inPos];
                bufferc=bufferc+bufferx*dfttC[i][j-inPos];
                
            }
            out->feat3[k/(4*GF_PASSFREQ/GF_PASSBND)][i]=sqrt(buffers*buffers+bufferc*bufferc);
        }
    }
    */

    
    // ********************************
    /*
    for (int k=0; k<GF_PASSFREQ; k++)
    {
        for (int l=inPos; l<inPos+m_intwin_total; l+=GF_PASSFREQ/GF_PASSBND)
        {
            bufferx = 0;
            for(int kk=l;kk<l+GF_PASSFREQ/GF_PASSBND;kk++) bufferx+=inBuf[l][kk];
            
            out->feat[k][(l-inPos)/(GF_PASSFREQ/GF_PASSBND)]=bufferx/(GF_PASSFREQ/GF_PASSBND);
        }
    }
     */
    // ********************************

    return true;
}

// *********************************************************************************************************

bool CFPExtractor::MakeFingerPrintColFromFV2(TFeatureVector2* inBuf, int inBufSize, int inBufPos, TFingerPrintCol2* out)
{
    if (inBufSize-inBufPos<2 || !inBuf || !out) return false;
    float tmp=0,tmp2=0,max=0,max1=0,max2=0,max3=0;
    
    int maxI=0,maxI1=0,maxI2=0,maxI3=0;
     
    for (int k=0; k<GF_PASSBND; k++)
    {
        tmp=tmp2=max=maxI=max1=maxI1=max2=maxI2=max3=maxI3=0;
        for(int i=3; i<GF_PASSFREQ2/2; i++)
        {
            if(max<inBuf[inBufPos].feat[k][i])
            {
                max=inBuf[inBufPos].feat[k][i];
                maxI=i;
            }
        }

        (*out)[k]=64+maxI;

    }

    float diff = 0, diff2 = 0, diff3 = 0, diff4 = 0;
    float norm1 = 0, norm2 = 0, norm3 = 0, cosdist =0;
    for (int k=0; k<GF_PASSBND; k++)
    {
        for(int i=0;i<GF_PASSFREQ2;i++)
        {
            diff+=(inBuf[inBufPos].feat[k][i]-inBuf[inBufPos+1].feat[k][i]);
            diff2+=(inBuf[inBufPos].feat[k][i]-inBuf[inBufPos+2].feat[k][i]);
            
            diff3+=(inBuf[inBufPos].feat[k][i]*inBuf[inBufPos+1].feat[k][i]);
            diff4+=(inBuf[inBufPos].feat[k][i]*inBuf[inBufPos+2].feat[k][i]);
            
            norm1+=powf(inBuf[inBufPos].feat[k][i],2);
            norm2+=powf(inBuf[inBufPos+1].feat[k][i],2);
            norm3+=powf(inBuf[inBufPos+2].feat[k][i],2);
        }
        
        if(diff>0.01)(*out)[k+GF_PASSBND]='U';
        else if(diff<-0.01)(*out)[k+GF_PASSBND]='D';
        else (*out)[k+GF_PASSBND]='E';
        diff=0;
        
        if(diff2>0.01)(*out)[k+2*GF_PASSBND]='U';
        else if(diff2<-0.01)(*out)[k+2*GF_PASSBND]='D';
        else (*out)[k+2*GF_PASSBND]='E';
        diff2=0;
        
        cosdist = diff3/(0.01+sqrt(norm1*norm2));
        if(cosdist>0&&cosdist<0.95)(*out)[k+3*GF_PASSBND]='A';
        if(cosdist>0.95&&cosdist<0.97)(*out)[k+3*GF_PASSBND]='B';
        if(cosdist>0.97&&cosdist<0.98)(*out)[k+3*GF_PASSBND]='C';
        if(cosdist>0.98&&cosdist<0.99)(*out)[k+3*GF_PASSBND]='D';
        if(cosdist>0.99&&cosdist<0.992)(*out)[k+3*GF_PASSBND]='E';
        if(cosdist>0.992&&cosdist<0.994)(*out)[k+3*GF_PASSBND]='F';
        if(cosdist>0.994&&cosdist<0.996)(*out)[k+3*GF_PASSBND]='G';
        if(cosdist>0.996&&cosdist<0.998)(*out)[k+3*GF_PASSBND]='H';
        if(cosdist>0.998&&cosdist<0.999)(*out)[k+3*GF_PASSBND]='I';
        if(cosdist>0.999&&cosdist<1.1)(*out)[k+3*GF_PASSBND]='X';
        diff3 = diff4 = norm1 = norm2 = norm3 = 0;
        
    }
  
    
    return true;
}

// *********************************************************************************************************

bool CFPExtractor::MakeFingerPrintColFromFV(TFeatureVector* inBuf, int inBufSize, int inBufPos, TFingerPrintCol* out)
{
    if (inBufSize-inBufPos<2 || !inBuf || !out) return false;
    
    float a,b,c,e,g,kl;
    float e1,g1;
    float e2,g2;
    a = b = c = e = g = 0.0f;
    e1 = g1 = 0.0f;
    e2 = g2 = 0.0f;
    for (int k=0; k<GF_PASSFREQ; k++)
    {
        float v12 = pow(inBuf[inBufPos].var[k],2);
        float v22 = pow(inBuf[inBufPos+1].var[k],2);
        float m12 = pow(inBuf[inBufPos].mean[k],2);
        float vd = inBuf[inBufPos].var[k] - inBuf[inBufPos+1].var[k];
        float md = inBuf[inBufPos].mean[k] - inBuf[inBufPos+1].mean[k];
        float tc = pow(md,2) / (0.00001+m12);
        float te = md / (0.00001+inBuf[inBufPos].mean[k]);
//        float td = pow(vd,2) / (0.00001+v12);
        float tg = vd / (0.00001+inBuf[inBufPos].var[k]);

        a += v12; b += v22; c += tc; e += te; g += tg; 
        if (k>=1 && k<GF_PASSFREQ/2)
        {
    	    e1 += te; g1 += tg; 
        } else if (k>=GF_PASSFREQ/2) {
    	    e2 += te; g2 += tg; 
        }
    }
    kl = a/(b+0.0000000001) + b/(a+0.0000000001) + c*(1/(a+0.0000000001) + 1/(b+0.0000000001));
    (*out)[0] = (kl>2.09f)?'A':((kl>2.01f)?'B':'C');
    (*out)[1] = (e>0.1f)?'A':((e<-0.1f)?'B':'C');
    e = inBuf[inBufPos].maxmean.value - inBuf[inBufPos+1].maxmean.value ;
    (*out)[2] = (e>0.05f)?'A':((e<-0.05f)?'B':'C');
    (*out)[3] = (g>0.1f)?'A':((g<-0.1f)?'B':'C');
    (*out)[4] = (e1>0.1f)?'A':((e1<-0.1f)?'B':'C');
    e1 = inBuf[inBufPos].maxmeanLow.value - inBuf[inBufPos+1].maxmeanLow.value ;
    (*out)[5] = (e1>0.05f)?'A':((e1<-0.05f)?'B':'C');
    (*out)[6] = (g1>0.001f)?'A':((g1<-0.001f)?'B':'C');
    (*out)[7] = (e2>0.02f)?'A':((e2<-0.02f)?'B':'C');
    e2 = inBuf[inBufPos].maxmeanHi.value - inBuf[inBufPos+1].maxmeanHi.value ;
    (*out)[8] = (e2>0.05f)?'A':((e2<-0.05f)?'B':'C');
    (*out)[9] = (g2>0.001f)?'A':((g2<-0.001f)?'B':'C');

    a = b = c = e = g = 0.0f;
    int where=9;
    for (int k=2; k<GF_PASSFREQ/2; k+=2)
    {
        where++;
        if (pow(inBuf[inBufPos].mean[k]+inBuf[inBufPos].mean[k+1],2) - pow(inBuf[inBufPos].mean[k+2]+inBuf[inBufPos].mean[k+2+1],2) - 
    	    ( pow(inBuf[inBufPos+1].mean[k]+inBuf[inBufPos+1].mean[k+1],2)-pow(inBuf[inBufPos+1].mean[k+2]+inBuf[inBufPos+1].mean[k+2+1],2) ) >0)
            (*out)[where]='U';
        else
            (*out)[where]='D';
    }
    return true;
}

// *********************************************************************************************************
bool CFPExtractor::MakeFFTFromData(short* data, int dataLen, int dataPos, double* out)
{
    if (!data || dataLen-dataPos<m_win_size || !out) return false;
    for (int i=0; i<m_win_size; i++)
	m_tmp_buf[i] = (double)data[i+dataPos]*m_hamming[i];
    m_FFT->make_fft(out, m_tmp_buf);
    int j = m_win_size/2;
    for (int i=0; i<j; i++)
    {
	out[i]/=(double)m_win_size;
	out[i+j]/=(double)m_win_size;
	if (fabs(out[i])>0.01 && fabs(out[i+j])>0.01)
	{
	    out[i+j] = (out[i]*out[i]+out[i+j]*out[i+j]);
	    out[i] = sqrt(out[i+j]);
	}
	else
	{
	    out[i] = 0.0;
	    out[i+j] = 0.0;
	}
    }
    m_sliding_norm_ptr++;
    if (m_sliding_norm_ptr>=GF_SLIDING_NORM_SIZE) { m_sliding_norm_ptr = 0; m_sliding_norm_over = true; }
    
    double norm = 0.0, norm1 = 0.0;
    for (int i=0; i<128; i++)
    {
	norm+=out[i+j];
	norm1+=pow(m_filt[i]*out[i],2);
    }
    m_sliding_norm[m_sliding_norm_ptr] = norm;
    m_sliding_norm1[m_sliding_norm_ptr] = norm1;

    return true;
}

// *********************************************************************************************************

bool CFPExtractor::MakeXFromFFT(double *fftbuf, TXVector* out)
{
    if (!fftbuf || !out) return false;
    double norm = 0.0, norm1 = 0.0;
    for (int i=0; i<GF_SLIDING_NORM_SIZE; i++)
    {
	norm += m_sliding_norm[i];
	norm1 += m_sliding_norm1[i];
    }
    double ratio = (norm1>0.0001)? norm/norm1 : 1.0;
//    printf("%f ",ratio);
    for (int i=0; i<128; i++)
	fftbuf[i] = ratio * fftbuf[i] * m_filt[i];
	
    for (int k=0; k<GF_PASSFREQ; k++)
    {
	float tmp = 0.0f;
	for (int kk = 3*k; kk<3*(k+1); kk++)
	    tmp+=fftbuf[kk];
	(*out)[k] = tmp;
    }
    return true;
}

// *********************************************************************************************************

int CFPExtractor::TakePossibleFV()
{
    if (m_cnt_featbuf<2) return 0;
    int np = m_cnt_featbuf-1;
    if (np > GF_WORKBUF_SIZE-m_cnt_resultbuf) np = GF_WORKBUF_SIZE-m_cnt_resultbuf;
    if (np<=0) return 0;

    for (int i=0; i<np; i++)
    {
        //MakeFingerPrintColFromFV(m_featbuf, m_cnt_featbuf, i, &(m_resultbuf[m_cnt_resultbuf++]));
        MakeFingerPrintColFromFV2(m_featbuf2, m_cnt_featbuf, i, &(m_resultbuf2[m_cnt_resultbuf++]));
    }
    return np;
}

// *********************************************************************************************************

int CFPExtractor::TakePossibleX()
{
    if (m_cnt_xbuf<m_intwin_total) return 0;
    int np = (m_cnt_xbuf-m_intwin_total+m_frac_used)/m_frac_used;
    if (np > GF_WORKBUF_SIZE-m_cnt_featbuf) np = GF_WORKBUF_SIZE-m_cnt_featbuf;
    if (np<=0) return 0;

    for (int i=0; i<np; i++)
    {
        //MakeFeatureVectorFromX(m_xbuf, m_cnt_xbuf, i*m_frac_used, &(m_featbuf[m_cnt_featbuf++]));
        MakeFeatureVector2FromX(m_xbuf, m_cnt_xbuf, i*m_frac_used, &(m_featbuf2[m_cnt_featbuf++]));
    }
    int a = TakePossibleFV();
    if (a>0 && a<=m_cnt_featbuf)
    {
	//memmove(&(m_featbuf[0]),&(m_featbuf[a]),sizeof(TFeatureVector)*(m_cnt_featbuf-a));
    memmove(&(m_featbuf2[0]),&(m_featbuf2[a]),sizeof(TFeatureVector2)*(m_cnt_featbuf-a));
	m_cnt_featbuf-=a;
    }
    return np*m_frac_used;
}

// *********************************************************************************************************

int CFPExtractor::TakePossibleData(short* data, int len)
{
    if (!data || len<m_win_size) return 0;

    int np = (len-m_win_size+m_win_step)/m_win_step;
    if (np > GF_WORKBUF_SIZE-m_cnt_xbuf) np = GF_WORKBUF_SIZE-m_cnt_xbuf;
    if (np<=0) return 0;
    
    if (!m_sliding_norm_over)
    {
	for (int i=0; i<np; i++)
	    MakeFFTFromData(data, len, i*m_win_step, m_fftbuf[i]);
	for (int i=0; i<np; i++)
	    MakeXFromFFT(m_fftbuf[i], &(m_xbuf[m_cnt_xbuf++]));
    } else {
	for (int i=0; i<np; i++)
	{
	    MakeFFTFromData(data, len, i*m_win_step, m_fftbuf[i]);
	    MakeXFromFFT(m_fftbuf[i], &(m_xbuf[m_cnt_xbuf++]));
	}	
    }
    int a = TakePossibleX();
    if (a>0 && a<=m_cnt_xbuf)
    {
	memmove(&(m_xbuf[0]),&(m_xbuf[a]),sizeof(TXVector)*(m_cnt_xbuf-a));
	m_cnt_xbuf-=a;
    }
    return np*m_win_step;
}

// *********************************************************************************************************

int CFPExtractor::GivePossibleResultSize()
{
    return m_cnt_resultbuf;
}

// *********************************************************************************************************
int CFPExtractor::GivePossibleResult(TFingerPrintCol* dest, int destSize)
{
    int np = m_cnt_resultbuf;
    if (np<destSize) np = destSize;
    if (np<=0 && !dest) return 0;
    
    memcpy(dest, m_resultbuf, sizeof(TFingerPrintCol)*np);
    if (np<m_cnt_resultbuf)
	memmove(&(m_resultbuf[0]), &(m_resultbuf[np]), sizeof(TFingerPrintCol)*(m_cnt_resultbuf-np));
    m_cnt_resultbuf-=np;

    return np;
}

// *********************************************************************************************************
int CFPExtractor::FillPossibleResult(char** fingerPrint, int destSize, int fromPos)
{
    int np = destSize-fromPos;
    if (np>m_cnt_resultbuf) np = m_cnt_resultbuf;
    if (np<=0 || !fingerPrint) return 0;
    for (int i=0; i<19; i++)
	if (!fingerPrint[i]) return 0;

    for (int i=0; i<np; i++)
	for (int j=0; j<19; j++)
	    fingerPrint[j][fromPos+i]=m_resultbuf[i][j];

    if (np<m_cnt_resultbuf)
	memmove(&(m_resultbuf[0]), &(m_resultbuf[np]), sizeof(TFingerPrintCol)*(m_cnt_resultbuf-np));
    m_cnt_resultbuf-=np;
    
    return np;
}

// *********************************************************************************************************
int CFPExtractor::FillPossibleResult2(char** fingerPrint, int destSize, int fromPos)
{
    int np = destSize-fromPos;
    if (np>m_cnt_resultbuf) np = m_cnt_resultbuf;
    if (np<=0 || !fingerPrint) return 0;
    for (int i=0; i<GF_FINELEM; i++)
        if (!fingerPrint[i]) return 0;
    
    for (int i=0; i<np; i++)
        for (int j=0; j<GF_FINELEM; j++)
            fingerPrint[j][fromPos+i]=m_resultbuf2[i][j];
    
    if (np<m_cnt_resultbuf)
        memmove(&(m_resultbuf2[0]), &(m_resultbuf2[np]), sizeof(TFingerPrintCol2)*(m_cnt_resultbuf-np));
    m_cnt_resultbuf-=np;
    
    return np;
}
