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


#ifndef GFINGERPRINT_H
#define GFINGERPRINT_H

#include "GFPExtractor.h"

#define CFP_OK			0
#define CFP_ERR_NOFILE		1
#define CFP_ERR_BADFORMAT	2
#define CFP_ERR_TOOSHORT	3

#define ENABLE_SIMD

extern int CF_SLOW_COMPARE_L;
extern int CF_SLOW_COMPARE_H;

class CFingerPrint
{
public:
    char *fingerprint[19];
    char *fingerprint2[GF_FINELEM];
    char name[255];
    int N;

    CFingerPrint();
    CFingerPrint(int n, char *fp[19]);
    //CFingerPrint(int n, char *fp[GF_PASSFREQ]);
    void fill(int n, char *fp[19]);
    void fill (int n, char *fp[19], int step);
    void fill2(int n, char *fp[GF_FINELEM]);
    void fill2 (int n, char *fp[GF_FINELEM], int step);
    bool save(const char *fname);
    bool save2(const char *fname);
    int open(const char *fname);
    int open2(const char *fname);
    int extract(const char* fname, CFPExtractor* extractor);
    int extract(const char* fname, CFPExtractor* extractor, int from, int length);
    int compare_to(CFingerPrint* dest, int &pos, bool full = true);
    int compare_to2(CFingerPrint* dest, int &pos, bool full = true);
    int compare_to_slow(CFingerPrint* dest, int &pos, bool full = true);
    ~CFingerPrint();
};

#endif
