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



#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <fcntl.h>
#include <string.h>
#include <vector>
#include <string>
#include <utility>
#include "GFingerPrint.h"



using namespace std;


int main(int argc, const char * argv[])
{
    if (argc<4)
    {
        printf("ERR\nMissing parameters\n");
        printf("\nUsage: gtag extract filename_audio filename_fingerprint (to extract fingerprint)\n");
        printf("\nUsage: gtag compare filename_audio directory_fingerprints (to compare filename to directory of fingerprints)\n");
        return 254;
    }
    
    if (strcmp(argv[1],"extract")==0)
    {
        // extracting fingerprint
        CFingerPrint *s = new CFingerPrint();
        int rez = s->extract(argv[2],NULL);
        if (rez!=CFP_OK)
        {
            printf("ERR\nCannot open file (code=%d)\n",rez);
            return 253;
        }
        if (!s->save2((char*)argv[3]))
        {
            printf("ERR\nCannot create file\n");
            return 251;
        }
        printf("OK\n");
        delete s;
        return 0;
    }
    
    if (strcmp(argv[1],"compare")==0)
    {
        int thold = 700;
        if (argc>4) thold = atoi(argv[4]);
        
        DIR *dp;
        struct dirent *dirp;
        if ((dp = opendir((char*)argv[3])) == NULL)
        {
            printf("ERR\nError reading fingerprint directory\n");
            return 252;
        }
        
        vector <string> list;
        
        while ((dirp = readdir(dp)) != NULL)
        {
            int l = (int)strlen(dirp->d_name);
            if (l>4)
            {
                if (strcmp(dirp->d_name+(l-4),".fin")==0)
                    list.push_back(string(dirp->d_name));
                
            }
        }
        closedir(dp);
        
        
        CFingerPrint *s = new CFingerPrint();
        CFingerPrint *s1 = new CFingerPrint();
        
        if (s->extract(argv[2],NULL)!=CFP_OK)
        {
            printf("ERR\nCannot open file\n");
            return 253;
        }
        
        int best_dist = 0;
        int best_i = -1;
        int best_pos = 0;
        int dist,disti;
        
        string fn;
        for (int i = 0; i<list.size(); i++) {
            fn = string(argv[3])+"/"+list[i];
            if (s1->open2(fn.c_str())!=CFP_OK)
                continue;
            dist = s->compare_to2(s1, disti, true);
            if (dist>best_dist)
            {
                best_dist = dist;
                best_i = i;
                best_pos = disti;
            }

        }
        
        if (best_i==-1)
        {
            printf("ERR\nNothing found\n");
        }
        else
        {
            if (best_dist < thold)
            {
                best_dist = 0;
            }
            if(best_dist > thold )
            {
                printf("OK\n%s\n%i,%i\n",list[best_i].c_str(),best_dist,best_pos);
            }
            
            else printf("OK\nNot recognized\n%i,%i\n",0,-1);
        }
        
        delete s;
        delete s1;
        return 0;
    }
    else
    {
        printf("ERR\nBad command parameter\n");
        return 254;
    }
}
