//
// Created by miken on 10/5/2018.
//
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include "nwops.h"

/* Docstrings */
static char version_str[] = "Version 0.1 (beta)";
static char module_docstring[] =
        "This module provides some Needleman-Wunch operations via C.";
static char print_string_lens_docstring[] =
        "Simply prints the int values provided for these strings, as a test";
static char get_version_docstring[] =
        "Just prints the version number of this module. For easy testing that import is ok.";
static char nw_most_matches_docstring[] =
        "Gets the NW-Score where m=1 and all other values are 0. Equivalent identifying the\n"
        "largest possible set of matches betwen the sequences";
static char get_best_startpos_docstring[] =
        "Uses the NW score (d=-1,m=1,g=-2) to decide which of the three sequences is best. \n"
        "                                                                                  \n"
        "Usage: find_best_start_pos(seq1, seq2, seq3, seq_protein) where seq1-seq3 are each\n"
        "of the 3 frame-shifted translations of the nucleotide sequence, and 'seq_protein' \n"
        "is the actual amino acid sequence. Returns a tuple of:                            \n"
        "   (frame_pos, match_ct, gap_ct, diff_ct, maxres, aln1, alnP)                     \n"
        "   where frame_pos is the frame-shift value (0, 1 or 2), max res is the score of  \n"
        "   best result, match_ct/gap_ct/diff_ct are the match/gap/difference counts of the\n"
        "   best result, and finally aln1 and alnP are the aligned versions of the         \n"
        "   sequences.";
static char get_pairwise_alignment_docstring[] =
        "retrieves the actual pairwise alignment for the two strings under a basic NW.     \n"
        "                                                                                  \n"
        "Usage:    get_pairwise_alignment(seq1, seq2, Fmat, Dmat)                          \n"
        "    where seq1, seq2 are strings that should be aligned. Fmat and Dmat should be  \n"
        "    2d Numpy arrays of type 'int32' and should be of size ( len(seq1) + 1,        \n"
        "    len(seq2) + 1). During the routine, Fmat is populated with the DP matrix      \n"
        "    computed during the NW algorithm. Dmat is populated with the \"arrows\" to get\n"
        "    the optimal path through the matrix, from top right to bottom left, where     \n"
        "    0 = from the left, 1 = from diagonal, and 2 = from above.                     \n"
        "                                                                                  \n"
        "Returns a tuple of:                                                               \n"
        "    (int) alignment length                                                        \n"
        "    (int) number of gaps                                                          \n"
        "    (int) best NW score                                                           \n"
        "    (str) aligned Sequence 1                                                      \n"
        "    (str) aligned Sequence 2                                                      \n"
        "";
static char test_numpy_data_exchange_docstring[] =
        "This just takes some data in the form of numpy arrays, makes some small changes, \n"
        "and then releases them, then checking to see if the changes persist. This was    \n"
        "a function I wrote for testing a while ago and didn't delete.                    \n"
        "\n"
        "Usage: test_numpy_data_exchange(int_arr, doub_arr), where int_arr is a 2-d numpy \n"
        "array of type int64, and doub_arr is a 2-d numpy array of type float64. This     \n"
        "function currently just adds 1 to every entry in the int_arr, and separately     \n"
        "prints the shapes of the two arrays.";

/* Available functions */
static PyObject *nwops_testPrintStringSizes(PyObject *self, PyObject *args);
static PyObject *nwops_getVersionString(PyObject *self, PyObject *args);
static PyObject *nwops_getMostMatches(PyObject *self, PyObject *args);
static PyObject *nwops_getStartFrameByMatchCt(PyObject *self, PyObject *args);
static PyObject *nwops_getPairwiseAlignment(PyObject *self, PyObject *args);
static PyObject *nwops_getPairwiseAlignment_testOpt(PyObject *self, PyObject *args);
static PyObject *nwops_testNumpyDataExchange(PyObject *self, PyObject *args);
static PyObject *nwops_change_data(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
        {"print_string_lens", nwops_testPrintStringSizes, METH_VARARGS, print_string_lens_docstring},
        {"get_version", nwops_getVersionString, METH_VARARGS, get_version_docstring},
        {"nw_most_matches", nwops_getMostMatches, METH_VARARGS, nw_most_matches_docstring},
        {"find_best_start_pos", nwops_getStartFrameByMatchCt, METH_VARARGS, get_best_startpos_docstring},
        {"get_pairwise_alignment", nwops_getPairwiseAlignment, METH_VARARGS, get_pairwise_alignment_docstring},
//        {"get_pairwise_alignment_test", nwops_getPairwiseAlignment_testOpt, METH_VARARGS, get_pairwise_alignment_numpy_docstring},
//        {"test_numpy_data_exchange", nwops_testNumpyDataExchange, METH_VARARGS, test_numpy_data_exchange_docstring},
        {"test_numpy_data_mod", nwops_change_data, METH_VARARGS, test_numpy_data_exchange_docstring},
//        {"pre_process_order_children", treeops_preProcOrderChildren, METH_VARARGS, pre_proc_order_children_docstring},
//        {"get_version", treeops_getVersionString, METH_VARARGS, get_version_docstring},
//        {"relocate_tree_by_deflection_or_edge_angle", treeops_relocateTreeByDeflectOrEdge, METH_VARARGS, relocate_by_deflection_edge_docstring},
//        {"right_turn_angle", treeops_rightTurnAngle, METH_VARARGS, right_turn_angle_docstring},
//        {"rotational_angle_to_ray", treeops_getRotationalAngleToRay, METH_VARARGS, rotational_angle_to_ray_docstring},
        {NULL, NULL, 0, NULL}
};

/* Initialize the module */
static struct PyModuleDef nwops =
{
    PyModuleDef_HEAD_INIT,
    "nwops", /* name of module */
    module_docstring, /* module documentation, may be NULL */
    -1,   /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
    module_methods
};

PyMODINIT_FUNC PyInit_nwops(void)
{


    PyObject *module = PyModule_Create(&nwops);

    /* Load `numpy` functionality. */
    import_array();
    return module;
}

static PyObject *nwops_getVersionString(PyObject *self, PyObject *args)
{
    PyObject *ret = Py_BuildValue("s", version_str);
    return ret;
}

static PyObject *nwops_testPrintStringSizes(PyObject *self, PyObject *args)
{
    char *s1;
    char *s2;
    int n1, n2;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "s#s#", &s1, &n1, &s2, &n2))
        return NULL;

    printf("n1=%d\tn2=%d\n",n1,n2);

    /* Build the output tuple */
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *nwops_getMostMatches(PyObject *self, PyObject *args)
{
    // Find the largest number of pairwise matches that can be found
    // between the two sequences..
    //
    long d, m, g, a, res;    // (d)ifference, (m)atch, (g)ap, (a)sterisk
    d = 0; m = 1; g = 0; a = 0;

    char *s1;
    char *s2;
    int n1, n2;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "s#s#", &s1, &n1, &s2, &n2))
        return NULL;

    res = needleScoreFromCosts(n1, s1, n2, s2, d, m, g, a);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("l", res);
    return ret;
}

static PyObject *nwops_getStartFrameByMatchCt(PyObject *self, PyObject *args)
{
    // Finds which of the three frame sequences has the best alignemtn
    // Returns now just the match count but the frame position that that inmplies.
    //
    long d, m, g, a;
//    long res1, res2, res3;
    long maxres;    // (d)ifference, (m)atch, (g)ap, (a)sterisk
    d = -1; m = 1; g = -2; a = 0;

    char *sprot;
    char *s1;
    char *s2;
    char *s3;
    int np, n1, n2, n3, j;
    long frame_pos, match_ct, gap_ct, diff_ct, maxlen, totlen;
    frame_pos = -1;
    match_ct = -1; gap_ct = -1; diff_ct = -1;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "s#s#s#s#", &s1, &n1, &s2, &n2, &s3, &n3, &sprot, &np))
        return NULL;
    totlen=maxOfThreeLongs(n1, n2, n3)+np;
    char* aln1 = (char *)malloc(totlen + 1);
    char* alnP = (char *)malloc(totlen + 1);
    for (j=0; j<totlen+1; j++) {
        aln1[j]=0;
        alnP[j]=0;
    }


    maxres = bestNeedleScoreOfThree(n1, s1, n2, s2, n3, s3, np, sprot, d, m, g, a, &frame_pos, &gap_ct, &match_ct,
                                    &diff_ct, aln1, alnP);


    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("lllllss", frame_pos, match_ct, gap_ct, diff_ct, maxres, aln1, alnP);
    free(aln1);
    free(alnP);
    return ret;
}

static PyObject *nwops_getPairwiseAlignment_old(PyObject *self, PyObject *args)
{
    // Computes the pairwise

    PyObject *fmat, *dmat;
    long d, m, g, a;
//    long maxres;    // (d)ifference, (m)atch, (g)ap, (a)sterisk
    d = 0; m = 1; g = 0; a = 0;

    char *s_prot;
    char *s_nuke; //also a protein, but from a nuke originally
    int l_prot, l_nuke;
    int n_gaps, aln_length;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "s#s#OO",&s_prot,&l_prot, &s_nuke, &l_nuke, &fmat, &dmat))
        return NULL;

    /* Interpret the input objects as numpy arrays. */
    PyObject *fmat_array = PyArray_FROM_OTF(fmat, NPY_LONG, NPY_ARRAY_INOUT_ARRAY);
    PyObject *dmat_array = PyArray_FROM_OTF(dmat, NPY_LONG, NPY_ARRAY_INOUT_ARRAY);

    /* If that didn't work, throw an exception. */
    if (fmat_array == NULL || dmat_array == NULL) {
        Py_XDECREF(dmat_array);
        Py_XDECREF(fmat_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    long *f_mat_ptr = (long*)PyArray_DATA(fmat_array);
    long *d_mat_ptr = (long*)PyArray_DATA(dmat_array);

    /* Call the external C function. */
    int totlen;
    int j = 0;
    long bestscore = 0;
    totlen= l_nuke + l_prot;
    char *aln_prot = (char *)malloc(totlen + 1);
    char *aln_nuke = (char *)malloc(totlen + 1);
    for (j=0; j<totlen+1; j++) {
        aln_prot[j]=0;
        aln_nuke[j]=0;
    }

//    long* myFmat = (long *)malloc((l_nuke+2) * (l_prot+2) * sizeof(long));
//    long* dirmat = (long *)malloc((l_nuke+2) * (l_prot+2) * sizeof(long));

    bestscore = getNeedleAlignment(l_prot, s_prot, l_nuke, s_nuke, f_mat_ptr, d, m, g, a,
                        aln_prot, aln_nuke, d_mat_ptr, &aln_length, &n_gaps);

    /* Clean up. */
    Py_DECREF(fmat_array);
    Py_DECREF(dmat_array);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("lllss", aln_length, n_gaps, bestscore, aln_prot, aln_nuke);


//    free(myFmat);
//    free(dirmat);
    free(aln_prot);
    free(aln_nuke);

    return ret;
}

//static PyObject *nwops_getPairwiseAlignment_toNumpy(PyObject *self, PyObject *args)
static PyObject *nwops_getPairwiseAlignment(PyObject *self, PyObject *args)
{
    // Computes the pairwise alignment of two input sequences. Optionally takes
    // two numpy array arguments that will be populated for analysis in Python.

    PyObject *fmat = NULL;
    PyObject *dmat = NULL;
    long d, m, g, a;
//    long maxres;    // (d)ifference, (m)atch, (g)ap, (a)sterisk
    d = 0; m = 1; g = 0; a = 0;

    char *s_prot;
    char *s_nuke; //also a protein, but from a nuke originally
    int l_prot, l_nuke;
    int n_gaps, aln_length;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "s#s#|OO",&s_prot,&l_prot, &s_nuke, &l_nuke, &fmat, &dmat))
        return NULL;

    long *f_mat_ptr;
    long *d_mat_ptr;

    PyObject *fmat_array;
    PyObject *dmat_array;

    if (fmat)
    {
        /* Interpret the input objects as numpy arrays. */
        fmat_array = PyArray_FROM_OTF(fmat, NPY_LONG, NPY_ARRAY_INOUT_ARRAY); //Interpret as numpy array
        if (fmat_array == NULL) {
            Py_XDECREF(fmat_array); // Kill and throw exception on failure.
            return NULL;
        }
        f_mat_ptr = (long*)PyArray_DATA(fmat_array); //Get pointer to data as C-type.
    } else {
        f_mat_ptr = (long *)malloc((l_nuke+2) * (l_prot+2) * sizeof(long));
    }

    if (dmat)
    {
        dmat_array = PyArray_FROM_OTF(dmat, NPY_LONG, NPY_ARRAY_INOUT_ARRAY);
        if (dmat_array == NULL) {
            Py_XDECREF(dmat_array);
            return NULL;
        }
        d_mat_ptr = (long*)PyArray_DATA(dmat_array);
    } else {
        d_mat_ptr = (long *)malloc((l_nuke+2) * (l_prot+2) * sizeof(long));
    }

    /* Call the external C function. */
    int totlen;
    int j = 0;
    long bestscore = 0;
    totlen= l_nuke + l_prot;
    char *aln_prot = (char *)malloc(totlen + 1);
    char *aln_nuke = (char *)malloc(totlen + 1);
    for (j=0; j<totlen+1; j++) {
        aln_prot[j]=0;
        aln_nuke[j]=0;
    }

    bestscore = getNeedleAlignment(l_prot, s_prot, l_nuke, s_nuke, f_mat_ptr, d, m, g, a,
                                   aln_prot, aln_nuke, d_mat_ptr, &aln_length, &n_gaps);

//    /* Clean up. */
    if (!fmat) {
        free(f_mat_ptr);
    } else {
        Py_XDECREF(fmat_array);
    }
    if (!dmat) {
        free(d_mat_ptr);
    } else {
        Py_XDECREF(dmat_array);
    }

//    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("lllss", aln_length, n_gaps, bestscore, aln_prot, aln_nuke);

    free(aln_prot);
    free(aln_nuke);

    return ret;
}

static PyObject *nwops_change_data(PyObject *self, PyObject *args)
{
    PyObject *int_mat_pyo, *doub_mat_pyo;


    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OO",&int_mat_pyo, &doub_mat_pyo))
    {
        return NULL;
    }
//    PyArrayObject  *int_mat_AO;
//    PyArrayObject  *doub_mat_AO;
    PyObject  *int_mat_AO;
    PyObject  *doub_mat_AO;

    /* Interpret the input objects as numpy arrays. */
    int_mat_AO = PyArray_FROM_OTF(int_mat_pyo, NPY_INT64, NPY_ARRAY_INOUT_ARRAY);
    doub_mat_AO = PyArray_FROM_OTF(doub_mat_pyo, NPY_FLOAT64, NPY_ARRAY_INOUT_ARRAY);

    /* If that didn't work, throw an exception. */
    if (int_mat_AO == NULL || doub_mat_AO == NULL ) {
        Py_XDECREF( int_mat_AO );
        Py_XDECREF( doub_mat_AO );
        return NULL;
    }

    /* How many data points are there? */
    int Nint1 = (int)PyArray_DIM(int_mat_AO, 0);
    int Nint2 = (int)PyArray_DIM(int_mat_AO, 1);
    int Ndoub1 = (int)PyArray_DIM(doub_mat_AO, 0);
    int Ndoub2 = (int)PyArray_DIM(doub_mat_AO, 1);
    printf("int matrix is: (%d , %d), doub matrix is (%d, %d)\n", Nint1, Nint2, Ndoub1, Ndoub2);

    /* Get pointers to the data as C-types. */
    long long *int_mat_ptr = (long long*)PyArray_DATA(int_mat_AO);
    double *doub_mat_ptr = (double*)PyArray_DATA(doub_mat_AO);

    int ct;
    ct = 0;
    int i, j;

    for (i=0; i<Nint1; i++) {
        for (j=0; j<Nint2; j++) {
            int_mat_ptr[i*Nint2+j]++;
        }
    }
    for (i=0; i<Ndoub1; i++) {
        for (j=0; j<Ndoub2; j++) {
            doub_mat_ptr[i*Ndoub2+j]+=0.5;
        }
    }

    Py_DECREF(doub_mat_AO);
    Py_DECREF(int_mat_AO);

    /* Build the output tuple */
    Py_INCREF(Py_None);
    return Py_None;
}
