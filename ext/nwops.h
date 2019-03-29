//
// Created by miken on 10/5/2018.
//

#ifndef EXT_NWOPS_H
#define EXT_NWOPS_H

long needleScore(int len1, char* seq1, int len2, char* seq2, long* F, long d, long m, long g, long a, long* dirmat);
long needleScoreFromCosts(int len1, char* seq1, int len2, char* seq2, long d, long m, long g, long a);
long bestNeedleScoreOfThree(int s1len, char* s1, int s2len, char* s2, int s3len,
                            char* s3, int s4len, char* s4, long d, long m, long g, long a, long *fr,
                            long *gap_ct, long *match_ct, long *diff_ct, char* aln1, char* aln4);
long which_maxOfThreeLongs(long res1, long res2, long res3);
long maxOfThreeLongs(long res1, long res2, long res3);
void intZeroOut(int *arr, long cells);
void longZeroOut(long *arr, long cells);
void inplace_reverse(char * str);
long getNeedleAlignment(int len1, char* seq1, int len2, char* seq2, long* F, long d, long m, long g,
                        long a, char* aln1, char* aln2, long* dirmat, int *alnLen, int *gapct);
void populate_position_array_from_string(char *aln_seq, int *pos_arr, int pos_arr_len);


#endif //EXT_NWOPS_H
