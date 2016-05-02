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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <dirent.h>
#include <fcntl.h>
#include <unistd.h>
#include "GFingerPrint.h"
#include "GFPExtractor.h"

int CF_SLOW_COMPARE_L=0;
int CF_SLOW_COMPARE_H=4;


static void fp_distance2 (int max1, int max2, char* query, char * data,int stepData, float & distOut,int &indOut,int &indOut1,int &indOut2)
{
    int *hit;
    
    if (max1<max2)
    {
        hit=new int[max2-max1+1];
        for (int i=0; i<max2-max1; i++) hit[i]=0;
        
        for (int k=0; k<max2-max1; k+=stepData)
        {
            for (int i=0; i<max1; i++)
            {
                if(std::abs(query[i]-data[i+k])<1)
                    hit[k]++;
            }
            //hit[k]=(max1*GF_PASSFREQ2-hit[k]);
        }
        
        int tmp,tmp1,tmp2;
        tmp=tmp1=tmp2=0;//hit[0];
        int ind,ind1,ind2;
        ind=ind1=ind2=0;
        for (int i=0; i<max2-max1; i++)
        {
            //printf("%i ",hit[i]);
            if (hit[i]>tmp) {
                tmp=hit[i];
                ind=i;
            }
        }
        for (int i=0; i<max2-max1; i++)
        {
            if ( (hit[i]>tmp1) && (hit[i]<tmp) ) {
                tmp1=hit[i];
                ind1=i;
            }
        }
        for (int i=0; i<max2-max1; i++)
        {
            if ( (hit[i]>tmp2) && (hit[i]<tmp1) ) {
                tmp2=hit[i];
                ind2=i;
            }
        }
        //	printf("NORM: %i:%i %i:%i %i:%i\n",(int)tmp,ind,(int)tmp1,(int)ind1,(int)tmp2,ind2);
        distOut= 40*(float)tmp/(float)max1;
        indOut=ind;
        indOut1=ind1;
        indOut2=ind2;
        delete []hit;
        
    }
    
    else {
        distOut=0;
        indOut=0;
        indOut1=0;
        indOut2=0;
    }
}

static void fp_distanceStatic2 (int offset,int max, char* query, char * data, float & distOut)
{
    int hit=0;
    //if(offset+max<time_used*frac_used)
    {
        for (int k=0; k<max; k++)
        {
            //hit+=fabs(query[k]-data[offset+k]);
            if (std::abs(query[k]-data[offset+k])<1)
                hit++;
        }
       // hit=(max*GF_PASSFREQ2-hit);
        
    }
    distOut= (float)hit/(float)max;
}

int Max (float *values,int offset)
{
    float maxx=values[0];
    int index=0;
    for (int i=0; i<offset; i++)
    {
        if (values[i]>maxx)
        {
            maxx=values[i];
            index=i;
        }
    }
    return index;
}

int Max (int *values,int offset)
{
    int maxx=values[0];
    int index=0;
    for (int i=0; i<offset; i++)
    {
        if (values[i]>maxx)
        {
            maxx=values[i];
            index=i;
        }
    }
    return index;
}


int fp_compare2(CFingerPrint *s1, CFingerPrint *s2, int &startat, bool HD)

{
    
    if (s1->N>s2->N || s1->N==0 || s2->N==0) return -1; // s2 is necessarly bigger than s1 since we search for s1 in s2
    
    int indO,indO1,indO2;
    float bufD[230],dista=0;
    int bufI[230];
    int where=0;
    int from = (HD)?0:CF_SLOW_COMPARE_L;
    int to = (HD)?GF_FINELEM:CF_SLOW_COMPARE_H;
    
    for (int i=from; i<to; i+=1)
    {
        float dist;
        fp_distance2(s1->N,s2->N,s1->fingerprint2[i],s2->fingerprint2[i],1,dist,indO,indO1,indO2);
        bufI[where]=indO;
        where++;
        bufI[where]=indO1;
        where++;
        bufI[where]=indO2;
        where++;
        dista+=dist ;
    }
    //dista/=(to-from);
    //printf("%f \n",dista);
    
    int lastwhere=where;
    where=0;
    dista=0;
    for (int R=0; R<lastwhere; R++)
    {
        for (int i=0; i<GF_FINELEM; i++)
        {
            float dist;
            fp_distanceStatic2 (bufI[R],s1->N,s1->fingerprint2[i],s2->fingerprint2[i],dist);
            //printf("%i/%i->%f\n",R,i,dist);
            dista+=dist/(GF_FINELEM-0);
        }
        bufD[where]=dista;
        dista=0;
        where++;
    }
    int L=Max (bufD,where);
    dista=int(bufD[L]*120*(16));
    startat=bufI[L];
    
    return (int)dista;
}



CFingerPrint::CFingerPrint()
{
    N=0;
    for (int i=0; i<GF_FINELEM; i++)
    {
        fingerprint2[i]=NULL;
    }
//    strcpy(name,"NONAME");
}


void CFingerPrint::fill2 (int n, char *fp[GF_FINELEM])
{
    for (int i=0; i<GF_FINELEM; i++)
    {
        if (fingerprint2[i]) {
            delete fingerprint2[i];
            fingerprint2[i] = NULL;
        }
    }
    N = n;
    for (int i=0; i<GF_FINELEM; i++)
    {
        fingerprint2[i] = new char[N+1];
        memcpy(fingerprint2[i],fp[i],n);
    }
}

void CFingerPrint::fill2 (int n, char *fp[GF_FINELEM], int step)
{
    for (int i=0; i<GF_FINELEM; i++)
    {
        if (fingerprint2[i]) {
            delete fingerprint2[i];
            fingerprint2[i] = NULL;
        }
    }
    N = n/step;
    for (int i=0; i<GF_FINELEM; i++)
    {
        fingerprint2[i] = new char[N+1];
        for (int j=0, p=0; j<N; j++, p+=step)
    	    fingerprint2[i][j]=fp[i][p];
    }
}


bool CFingerPrint::save2(const char *fname)
{
    FILE* f = fopen (fname,"wb");
    if (!f) return false;
    fwrite(&N, 1, 4, f);
   // char ch;
    char ch = '\n';
    for (int s=0; s<GF_FINELEM; s++)
    {
        fwrite (this->fingerprint2[s],1,N,f);
        fwrite (&ch, 1,1, f);
    }
    fclose (f);
    return true;
}



int CFingerPrint::open2(const char *fname)
{
    struct flock fl;
    
    fl.l_type   = F_RDLCK;  /* F_RDLCK, F_WRLCK, F_UNLCK    */
    fl.l_whence = SEEK_SET;
    fl.l_start  = 0;
    fl.l_len    = 0;
    fl.l_pid    = getpid();
    FILE* f = fopen (fname,"rb");
    if (!f) return CFP_ERR_NOFILE;
    fcntl(fileno(f), F_SETLKW, &fl);
    
    fseeko(f, 0, SEEK_END);
    int fsize= (int)ftello (f);
    fseeko (f, 0, 0);
    
    char* buf = new char[fsize];
    fread ( (void *)buf,1,fsize,f);
    fclose (f);
    
    
    int n = ((int*)buf)[0];
    if (n<=0 || n>20000000) {
        delete[] buf;
        return CFP_ERR_BADFORMAT;
    }
    N = n;
    for (int i=0; i<GF_FINELEM; i++)
    {
        if (fingerprint2[i]) delete fingerprint2[i];
        fingerprint2[i] = new char[N+1];
        memcpy(fingerprint2[i],buf+4+i*(N+1),N);
    }
    delete[] buf;
    return CFP_OK;
}


int CFingerPrint::extract(const char* fname, CFPExtractor* extractor)
{
    return extract(fname, extractor, 0, 0);
}


int CFingerPrint::extract(const char* fname, CFPExtractor* extractor, int from, int length)
{
    FILE* f = fopen(fname,"rb");
    if (!f) return CFP_ERR_NOFILE;
    unsigned int fs1=0,ds1=0,sr1=0,sc1=0;
    fseek(f,4,0);
    fread(&fs1,1,4,f);
    fseek(f,0x14,0);
    fread(&sc1,1,4,f);
    fread(&sr1,1,4,f);
    int i;
    int dp1=0;
    for (i=0x20; i<0x4000; i++)
    {
	fseek(f,i,0);
	fread(&ds1,1,4,f);
	if (ds1=='atad') {
	    dp1=i;
	    break;
	}
    }
    if (dp1==0) {
	fclose(f);
	return CFP_ERR_BADFORMAT;
    }
    if (sr1!=8000 || sc1!=65537) {
	fclose(f);
	return CFP_ERR_BADFORMAT;
    }
    int l = 0;
    fread(&l,1,4,f);
    if (l==0) {
	// workaround
	fseek(f,0,2);
	l = (int)ftell(f)-dp1-8;
	fseek(f,dp1+8,0);
    }
    if (l<16000) {
	fclose(f);
	return CFP_ERR_TOOSHORT;
    }
    if (l>1600000000) l = 1600000000;
    l/=2;
    from*=8; // was in milliseconds
    length*=8;
    if (from>0 && from<l)
    {
	l-=from;
	fseek(f,dp1+8+from*2,0);
//	printf("Starting from %d\n",dp1+8+from*2);
    }
    if (length>0 && length<l)
	l=length;
    
    int n = l/800; // for 8 kHz and 10 vectors per second
    
    for (i=0; i<GF_FINELEM; i++)
    {
        if (fingerprint2[i]) delete fingerprint2[i];
        fingerprint2[i] = new char[n];
    }
    CFPExtractor* local_extractor= NULL;
    if (extractor==NULL)
    {
	extractor = local_extractor = new CFPExtractor(256, 80, 100, 10, 0.2f);
    }
    int BS = 100000;
    short* buf = new short[BS];
    int lptr = 0;
    int fptr = 0;
    int t = 0;
    while (l>0 || t>0)
    {
	int r = (l>(BS-lptr))?(BS-lptr):l;
	l-=r;
	fread(buf+lptr,2,r,f);
	t = extractor->TakePossibleData(buf, lptr+r);
	//int e = extractor->FillPossibleResult(fingerprint, n, fptr);
    int e = extractor->FillPossibleResult2(fingerprint2, n, fptr);
	fptr += e;
	memmove(buf, buf+t, (lptr+r-t)*2);
	lptr = (lptr+r-t);
    }
    if (local_extractor) delete local_extractor;
    delete[] buf;
    N = fptr;
    return CFP_OK;
}


int CFingerPrint::compare_to2(CFingerPrint* dest, int &pos, bool full)
{
    if (!dest || N==0) return -1;
    if (dest->N==0) return -1;
    int r = fp_compare2(this, dest, pos, full);
    pos = pos*100;
    return r;
}



CFingerPrint::~CFingerPrint()
{
    for (int i=0; i<GF_FINELEM; i++)
    {
        if (fingerprint2[i]) delete fingerprint2[i];
    }
}
