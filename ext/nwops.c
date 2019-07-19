//
// Created by miken on 10/5/2018.
//
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "nwops.h"

void longZeroOut(long *arr, long cells)
{
    int i;
    for (i=0;i<cells; i++)
    {
        arr[i]=0;
    }
}
void intZeroOut(int *arr, long cells)
{
    int i;
    for (i=0;i<cells; i++)
    {
        arr[i]=0;
    }
}

long maxOfThreeLongs(long res1, long res2, long res3)
{
    long maxres;
    if (res1 > res2) {
        if (res1 > res3) {
            maxres = res1;
//            frame = 0;
        } else {
            maxres = res3;
//            frame = 2;
        }
    } else {
        if (res2 > res3) {
            maxres = res2;
//            frame = 1;
        } else {
            maxres = res3;
//            frame = 2;
        }
    }
    return maxres;
}

long which_maxOfThreeLongs(long res1, long res2, long res3)
{
    long frame;
    if (res1 > res2) {
        if (res1 > res3) {
            frame = 0;
        } else {
            frame = 2;
        }
    } else {
        if (res2 > res3) {
            frame = 1;
        } else {
            frame = 2;
        }
    }
    return frame;
}
/*
 *  needleScore:
 *      len1/seq1 - sequence 1 and it's length
 *      len2/seq2 - ...
 *      F         - Big DP matrix for scoring (len1 + 1) x (len2 + 1)
 *      d/m/g/a   - Costs for deletion/match/gap/asterix
 *      dirmat    - Matrix of directions from which the F-values were found
 *
 * */
long needleScore(int len1, char* seq1, int len2, char* seq2, long* F, long d, long m, long g, long a, long* dirmat)
{
    // assumes that Fmat is big enough to handle the scores
    // d = difference score, m = match score, g = gap, a = asterix
    int i, j;
    long mat, del, ins, tmp;

    int dir, dirtmp;
    // char *ast = "*";

    // fill the edges of F:
    for (i=0; i<len1+1;i++)
    {
        F[i*(len2+1)]=i*g; //should this be i*(len2+1)+1?
    }
    for (j=0; j<len2+1;j++)
    {
        F[j]=j*g;
    }
    // fill the middle of F:
    for (i=1; i<len1+1; i++)
    {
        for (j=1; j<len2+1; j++)
        {
            mat = F[(i-1)*(len2+1)+(j-1)] + (seq1[i-1]==42 || seq2[j-1]==42 ? a : (seq1[i-1]==seq2[j-1] ? m : d)); // (dir = 1)
            del = F[(i-1)*(len2+1)+(j)] + g; // from above (dir = 2)
            ins = F[(i)*(len2+1)+(j-1)] + g; // from the left (dir = 0)
            tmp = (mat > del ? mat : del);
            dirtmp = (mat > del ? 1 : 2);
            F[i*(len2+1)+(j)] = (tmp > ins ? tmp : ins);
            dir = (tmp > ins ? dirtmp : 0);
            dirmat[i*(len2+1)+(j)] = dir;
        }
    }

    return F[(len2+1)*(len1+1)-1];
}

long needleScoreFromCosts(int len1, char* seq1, int len2, char* seq2, long d, long m, long g, long a)
{
    long res; // (d)ifference, (m)atch, (g)ap, (a)sterisk
    long* myFmat = (long *)malloc((len1+2) * (len2+2) * sizeof(long));
    long* dirmat = (long *)malloc((len1+2) * (len2+2) * sizeof(long));

//    d = -1; m = 2; g = -4; a = 1;
//    d = 0; m = 1; g = 0; a = 0;

    res = needleScore(len1, seq1, len2, seq2, myFmat, d, m, g, a, dirmat);
    free(myFmat);
    free(dirmat);

    return(res);
}

long bestNeedleScoreOfThree(int s1len, char* s1, int s2len, char* s2, int s3len,
                            char* s3, int s4len, char* s4, long d, long m, long g, long a, long *fr,
                            long *gap_ct, long *match_ct, long *diff_ct, char* aln1, char* aln4)
{
    // (d)ifference, (m)atch, (g)ap, (a)sterisk
    long maxlen, res1, res2, res3, totlen;
    maxlen = maxOfThreeLongs(s1len, s2len, s3len);

    totlen = s4len + maxlen;
//    char* aln1 = (char *)malloc(totlen + 1);
    char* aln2 = (char *)malloc(totlen + 1);
    char* aln3 = (char *)malloc(totlen + 1);
//    char* aln4 = (char *)malloc(totlen + 1);
    char* aln5 = (char *)malloc(totlen + 1);
    char* aln6 = (char *)malloc(totlen + 1);
    int j;
    for (j=0; j<totlen+1; j++) {
        aln2[j]=0;
        aln3[j]=0;
        aln5[j]=0;
        aln6[j]=0;
    }

    long* myFmat = (long *)malloc((maxlen+2) * (s4len+2) * sizeof(long));
    long* dirmat = (long *)malloc((maxlen+2) * (s4len+2) * sizeof(long));

    int n_gaps, n_gaps2, n_gaps3, aln_length, aln_length2, aln_length3, maxres;
    match_ct[0] = 0;  diff_ct[0] = 0; gap_ct[0] = 0;

    res1 = getNeedleAlignment(s1len, s1, s4len, s4, myFmat, d, m, g, a, aln1, aln4, dirmat, &aln_length, &n_gaps);
    longZeroOut(myFmat, (maxlen+2) * (s4len+2));
    res2 = getNeedleAlignment(s2len, s2, s4len, s4, myFmat, d, m, g, a, aln2, aln5, dirmat, &aln_length2, &n_gaps2);
    longZeroOut(myFmat, (maxlen+2) * (s4len+2));
    res3 = getNeedleAlignment(s3len, s3, s4len, s4, myFmat, d, m, g, a, aln3, aln6, dirmat, &aln_length3, &n_gaps3);

    maxres=maxOfThreeLongs(res1, res2, res3);

    if (res2 > res1) {
        if (res3 > res2) {
            strcpy(aln1, aln3);
            strcpy(aln4, aln6);
            aln_length = aln_length3;
            n_gaps = n_gaps3;
            res1 = res3;
        } else {
            strcpy(aln1,aln2);
            strcpy(aln4,aln5);
            aln_length = aln_length2;
            n_gaps = n_gaps2;
            res1 = res2;
        }
    }

    int i;
    for (i=0; i<aln_length; i++) {
        if (aln1[i]==45 || aln4[i]==45) {
            gap_ct[0]++;
        } else if (aln1[i]==aln4[i]) {
            match_ct[0]++;
        } else {
            diff_ct[0]++;
        }
    }
    fr[0] = which_maxOfThreeLongs(res1, res2, res3);




    free(aln2);
    free(aln3);
    free(aln5);
    free(aln6);
    free(myFmat);
    free(dirmat);
    return(maxres);
}

//void getNeedleAlignemnt_Provide_F_dir(int len1, char* seq1, int len2, char* seq2, long* F, long d, long m, long g,
//                                      long a, char* aln1, char* aln2, int* dirmat, int *alnLen, int *gapct)

long getNeedleAlignment(int len1, char* seq1, int len2, char* seq2, long* F, long d, long m, long g,
                             long a, char* aln1, char* aln2, long* dirmat, int *aln_len, int *gapct)
{
    int p, r, c, n_gaps;
//    long mat, del, ins, tmp;
    long best_score;
    int dir;
//    int dirtmp;

    best_score = needleScore(len1, seq1, len2, seq2, F, d, m, g, a, dirmat);

    p = 0; n_gaps=0;
    r = len1; c = len2;
    // get the sequences
    while (r > 0 || c > 0)
    {
        dir = dirmat[r*(len2+1)+c];
        if (r==0) dir = 0;
        if (c==0) dir = 2;
        if (dir == 0)
        {
            aln1[p] = 45;
            aln2[p] = seq2[c-1];
            c--;
            n_gaps++;
        } else if (dir == 2) {
            aln1[p] = seq1[r-1];
            aln2[p] = 45;
            r--;
            n_gaps++;
        } else {
            aln1[p] = seq1[r-1];
            aln2[p] = seq2[c-1];
            r--;
            c--;
        }
        p++;
    }
    inplace_reverse(aln1);
    inplace_reverse(aln2);
    aln1[p] = 0;
    aln2[p] = 0;
    aln_len[0] = p;
    gapct[0]=n_gaps;
    return best_score;
}
void inplace_reverse(char * str)
{
    if (str)
    {
        char * end = str + strlen(str) - 1;

        // swap the values in the two given variables
        // XXX: fails when a and b refer to same memory location
#   define XOR_SWAP(a,b) do\
    {\
      a ^= b;\
      b ^= a;\
      a ^= b;\
    } while (0)

        // walk inwards from both ends of the string,
        // swapping until we get to the middle
        while (str < end)
        {
            XOR_SWAP(*str, *end);
            str++;
            end--;
        }
#   undef XOR_SWAP
    }
}

void populate_position_array_from_string(char *aln_seq, int *pos_arr, int pos_arr_len)
{
    // This is just going to run along the string and write down the position
    // whenever the character in the string is a non-gap. So for the string
    // 'ACT----GG' it will populate pos_arr to be: [0, 1, 2, -1, -1, -1, -1, 3, 4]
    //      aln_seq: string to be analyzed.
    //      pos_arr: array of positions, will be modified
    //      pos_arr_len: length of the variable pos_arr, to avoid a segfault.
    int loc_in_seq, loc_in_posarr,  aln_length;
    loc_in_posarr=0;
    loc_in_seq = 0;
    aln_length = strlen(aln_seq);

    for (loc_in_seq=0; loc_in_seq<aln_length; loc_in_seq++)
    {
        if (aln_seq[loc_in_seq]!=45)
        {
            if (loc_in_posarr>pos_arr_len) {
                printf("for some reason loc_in_posarr has gotten larger than pas_arr_len. (%d vs. %d, respectively",
                       loc_in_posarr, pos_arr_len);
                return;
            }
            pos_arr[loc_in_posarr]=loc_in_seq;
            loc_in_posarr++;
        }
        else
        {
            pos_arr[loc_in_posarr]=-1;
        }
    }
}