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



#ifndef GFPEXTRACTOR_H
#define GFPEXTRACTOR_H

#include "GFFT.h"

#define GF_WORKBUF_SIZE		2048
#define GF_PASSFREQ		40
#define GF_PASSBND		8
#define GF_FINELEM  32
#define GF_PASSFREQ2		40
#define GF_SLIDING_NORM_SIZE	200

typedef struct
{
    float	value;
    int		index;
} TPitch;

typedef float TDCTVector[128];
typedef float TXVector[GF_PASSFREQ];

typedef struct {
    float	mean[GF_PASSFREQ];
    float	var[GF_PASSFREQ];
    TPitch	maxmeanLow;
    TPitch	maxmeanHi;
    TPitch	maxmean;
} TFeatureVector;

typedef struct {
    double	feat[GF_PASSBND][GF_PASSFREQ2];
} TFeatureVector2;

typedef struct {
    short int	mat[GF_PASSFREQ/GF_PASSBND][GF_PASSFREQ/GF_PASSBND];
} TGabor;

typedef char TFingerPrintCol[19];
typedef char TFingerPrintCol2[GF_FINELEM];

class CFPExtractor
{
private:
    GFFT	*m_FFT;
    float	*m_ITWhamming;
    float *dfttS[GF_PASSFREQ2];
    float *dfttC[GF_PASSFREQ2];
    TXVector	m_xbuf[GF_WORKBUF_SIZE];
    TFeatureVector m_featbuf[GF_WORKBUF_SIZE];
    TFeatureVector2 m_featbuf2[GF_WORKBUF_SIZE];
    TGabor Gbr[8];
    float Gr[GF_PASSFREQ][100];
    float Gi[GF_PASSFREQ][100];
    double	*m_fftbuf[GF_WORKBUF_SIZE];
    int		m_cnt_xbuf;
    int		m_cnt_featbuf;
    int		m_frac_used;
    int		m_intwin_size;
    int		m_win_size;
    int		m_win_step;
    int		m_intwin_total;
    int		m_intwin_overlap;
    TFingerPrintCol m_resultbuf[GF_WORKBUF_SIZE];
    TFingerPrintCol2 m_resultbuf2[GF_WORKBUF_SIZE];
    int		m_cnt_resultbuf;
    double	*m_hamming;
    double	*m_tmp_buf;
    double	*m_sliding_norm;
    double	*m_sliding_norm1;
    int 	m_sliding_norm_ptr;
    bool	m_sliding_norm_over;
    double	m_filt[256];
// functions
    int TakePossibleX();
    int TakePossibleFV();
    bool MakeFFTFromData(short* data, int dataLen, int dataPos, double* out);
    bool MakeXFromFFT(double *fftbuf, TXVector* out);
    bool MakeFingerPrintColFromFV(TFeatureVector* inBuf, int inBufSize, int inBufPos, TFingerPrintCol* out);
    bool MakeFingerPrintColFromFV2(TFeatureVector2* inBuf, int inBufSize, int inBufPos, TFingerPrintCol2* out);
    bool MakeFeatureVectorFromX(TXVector* inBuf, int inBufSize, int inPos, TFeatureVector* out);
    bool MakeFeatureVector2FromX(TXVector* inBuf, int inBufSize, int inPos, TFeatureVector2* out);

public:
    CFPExtractor(int inWinSize, int inWinStep, int inIntWinSize, int inFracUsed, float inOverlap);
    ~CFPExtractor();
    
    int TakePossibleData(short* data, int len);
    int GivePossibleResultSize();
    int GivePossibleResult(TFingerPrintCol* dest, int destSize);
    int FillPossibleResult(char** fingerPrint, int destSize, int fromPos);
    int FillPossibleResult2(char** fingerPrint, int destSize, int fromPos);
};
#endif
