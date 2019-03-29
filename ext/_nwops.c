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
        "Gets the NW-Score where m=1 and all other values are 0. Equivalent identifying the"
        "largest possible set of matches betwen the sequences";
static char get_best_startpos_docstring[] =
        "Uses the NW score (d=-1,m=1,g=-2) to decide which of the three sequences is best. "
        "Usage: find_best_start_pos(seq1, seq2, seq3, seq_protein) where seq1-seq3 are each"
        "of the 3 frame-shifted translations of the nucleotide sequence, and 'seq_protein' "
        "is the actual amino acid sequence. Returns a tuple of:"
        "   (frame_pos, match_ct, gap_ct, diff_ct, maxres, aln1, alnP)"
        "   where frame_pos is the frame-shift value (0, 1 or 2), max res is the score of  "
        "   best result, match_ct/gap_ct/diff_ct are the match/gap/difference counts of the"
        "   best result, and finally aln1 and alnP are the aligned versions of the "
        "   sequences.";
static char get_pairwise_alignment_docstring[] =
        "retrieves the actual pairwise alignment for the two strings under a basic NW";
static char get_pairwise_alignment_numpy_docstring[] =
        "retrieves the actual pairwise alignment for the two strings under a basic NW, but"
        "this version depends on getting arrays from numpy so the results can be analyzed"
        "in python later.";
static char test_numpy_data_exchange_docstring[] =
        "This just takes some data in the form of numpy arrays, makes some small changes,"
        "and then releases them. We are checking to see if the changes persist";

/* Available functions */
static PyObject *nwops_testPrintStringSizes(PyObject *self, PyObject *args);
static PyObject *nwops_getVersionString(PyObject *self, PyObject *args);
static PyObject *nwops_getMostMatches(PyObject *self, PyObject *args);
static PyObject *nwops_getStartFrameByMatchCt(PyObject *self, PyObject *args);
static PyObject *nwops_getPairwiseAlignment(PyObject *self, PyObject *args);
static PyObject *nwops_getPairwiseAlignment_toNumpy(PyObject *self, PyObject *args);
static PyObject *nwops_testNumpyDataExchange(PyObject *self, PyObject *args);

/* Module specification */
static PyMethodDef module_methods[] = {
        {"print_string_lens", nwops_testPrintStringSizes, METH_VARARGS, print_string_lens_docstring},
        {"get_version", nwops_getVersionString, METH_VARARGS, get_version_docstring},
        {"nw_most_matches", nwops_getMostMatches, METH_VARARGS, nw_most_matches_docstring},
        {"find_best_start_pos", nwops_getStartFrameByMatchCt, METH_VARARGS, get_best_startpos_docstring},
        {"get_pairwise_alignment", nwops_getPairwiseAlignment, METH_VARARGS, get_pairwise_alignment_docstring},
        {"get_pairwise_alignment_tonumpy", nwops_getPairwiseAlignment_toNumpy, METH_VARARGS, get_pairwise_alignment_numpy_docstring},
        {"test_numpy_data_exchange", nwops_testNumpyDataExchange, METH_VARARGS, test_numpy_data_exchange_docstring},
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


static PyObject *nwops_getPairwiseAlignment(PyObject *self, PyObject *args)
{
    // Finds which of the three frame sequences has the best alignemtn
    // Returns now just the match count but the frame position that that inmplies.

//    PyObject *p_positions, *n_positions;
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
    PyObject *fmat_array = PyArray_FROM_OTF(fmat, NPY_LONG, NPY_ARRAY_IN_ARRAY);
    PyObject *dmat_array = PyArray_FROM_OTF(dmat, NPY_LONG, NPY_ARRAY_IN_ARRAY);

    /* If that didn't work, throw an exception. */
//    if (p_pos_array == NULL || n_pos_array == NULL) {
    if (fmat_array == NULL || dmat_array == NULL) {
        Py_XDECREF(dmat_array);
        Py_XDECREF(fmat_array);
        return NULL;
    }

    /* How many data points are there? */
//    int p_arr_len = (int)PyArray_DIM(p_pos_array, 0);
//    int n_arr_len = (int)PyArray_DIM(n_pos_array, 0);

    /* Get pointers to the data as C-types. */
    long *f_mat_ptr = (long*)PyArray_DATA(fmat_array);
    long *d_mat_ptr = (long*)PyArray_DATA(dmat_array);

    /* Call the external C function to compute the new angles. */
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

static PyObject *nwops_getPairwiseAlignment_toNumpy(PyObject *self, PyObject *args)
{
    // This is going to execult the pairwise alignment for two sequences, but it
    // is going to use everything taken carefully from python so that it will return
    // the data for analysis.

    PyObject *p_str_pyo, *n_str_pyo, *Fmat_pyo, *dirmat_pyo;
//    PyObject *Fmat_pyo, *dirmat_pyo;
    long d, m, g, a;
    // (d)ifference, (m)atch, (g)ap, (a)sterisk
     d = 0; m = 1; g = 0; a = 0;

    char *s_prot;
    char *s_nuke; //also a protein, but from a nuke originally
    int l_prot, l_nuke;
    int n_gaps, aln_length;

    /* Parse the input tuple */
//    if (!PyArg_ParseTuple(args, "s#s#OOOOiiii",&s_prot,&l_prot, &s_nuke, &l_nuke, &p_str_pyo, &n_str_pyo,
//                            &Fmat_pyo, &dirmat_pyo, &d, &m, &g, &a))
    if (!PyArg_ParseTuple(args, "s#s#OOOO", &s_prot, &l_prot, &s_nuke, &l_nuke, &p_str_pyo, &n_str_pyo, &Fmat_pyo, &dirmat_pyo))
        return NULL;
//    if (!PyArg_ParseTuple(args, "s#s#OO",&s_prot,&l_prot, &s_nuke, &l_nuke, &Fmat_pyo, &dirmat_pyo))
//        return NULL;

    /* Interpret the input objects as numpy arrays. */
    PyArrayObject  *p_str_uint8 = PyArray_FROM_OTF(p_str_pyo, NPY_UINT8, NPY_ARRAY_IN_ARRAY);
    PyArrayObject  *n_str_uint8 = PyArray_FROM_OTF(n_str_pyo, NPY_UINT8, NPY_ARRAY_IN_ARRAY);
    PyArrayObject  *Fmat_AO = PyArray_FROM_OTF(Fmat_pyo, NPY_INT64, NPY_ARRAY_IN_ARRAY);
    PyArrayObject  *dirmat_AO = PyArray_FROM_OTF(dirmat_pyo, NPY_INT64, NPY_ARRAY_IN_ARRAY);

    /* If that didn't work, throw an exception. */
    if (p_str_uint8 == NULL || n_str_uint8 == NULL || Fmat_AO == NULL || dirmat_AO == NULL ) {
//    if (Fmat_AO == NULL || dirmat_AO == NULL ) {
        Py_XDECREF(p_str_uint8);
        Py_XDECREF(n_str_uint8);
        Py_XDECREF( Fmat_AO );
        Py_XDECREF( dirmat_AO );
        return NULL;
    }


    /* Get pointers to the data as C-types. */
        // int *p_str_ptr = (int*)PyArray_DATA(p_pos_array);
        // int *n_str_ptr = (int*)PyArray_DATA(n_pos_array);
    char *p_str_ptr = (char*)PyArray_DATA(p_str_uint8);
    char *n_str_ptr = (char*)PyArray_DATA(n_str_uint8);
    long *Fmat_ptr = (long*)PyArray_DATA(Fmat_AO);
    long *dirmat_ptr = (long*)PyArray_DATA(dirmat_AO);

    /* Call the external C function to compute the new angles. */
    int totlen;
    totlen= l_nuke + l_prot;
    char *aln_prot = (char *)malloc(totlen + 1);
    char *aln_nuke = (char *)malloc(totlen + 1);

//    long* myFmat = (long *)malloc((l_nuke+2) * (l_prot+2) * sizeof(long));
//    long* dirmat = (long *)malloc((l_nuke+2) * (l_prot+2) * sizeof(long));

    getNeedleAlignment(l_prot, s_prot, l_nuke, s_nuke, Fmat_ptr, d, m, g, a,
                       aln_prot, aln_nuke, dirmat_ptr, &aln_length, &n_gaps);

//    populate_position_array_from_string(aln_prot, p_position_ptr, p_arr_len);
//    populate_position_array_from_string(aln_nuke, n_position_ptr, n_arr_len);

    int i;
    for (i=0; i<aln_length; i++)
    {
        p_str_ptr[i] = aln_prot[i];
        n_str_ptr[i] = aln_nuke[i];
    }


//    free(myFmat);
//    free(dirmat);
    free(aln_prot);
    free(aln_nuke);

    /* Clean up. */
    Py_XDECREF(p_str_uint8);
    Py_XDECREF(n_str_uint8);
    Py_DECREF(dirmat_AO);
    Py_DECREF(Fmat_AO);

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("ll", aln_length, n_gaps);
    return ret;
}

static PyObject *chi2_move_data(PyObject *self, PyObject *args)
{
    PyObject *Fmat_pyo, *dirmat_pyo;
    printf("ok, we are starting the function at the top\n");

    if (!PyArg_ParseTuple(args, "OO",&Fmat_pyo, &dirmat_pyo))
    {
        return NULL;
    }
    PyArrayObject  *Fmat_AO;
    PyArrayObject  *dirmat_AO;
    printf("state 2\n");

    /* Interpret the input objects as numpy arrays. */
    Fmat_AO = PyArray_FROM_OTF(Fmat_pyo, NPY_INT64, NPY_ARRAY_IN_ARRAY);
    dirmat_AO = PyArray_FROM_OTF(dirmat_pyo, NPY_INT64, NPY_ARRAY_IN_ARRAY);

    /* If that didn't work, throw an exception. */
    if (Fmat_AO == NULL || dirmat_AO == NULL ) {
        Py_XDECREF( Fmat_AO );
        Py_XDECREF( dirmat_AO );
        return NULL;
    }
    printf("state 3\n");

     /* How many data points are there? */
    int Nr = (int)PyArray_DIM(Fmat_AO, 0);
    int Nc = (int)PyArray_DIM(Fmat_AO, 1);
    printf("Size of F is: (%d , %d)\n", Nr, Nc);

    /* Get pointers to the data as C-types. */
    long *Fmat_ptr = (long*)PyArray_DATA(Fmat_AO);
    long *dirmat_ptr = (long*)PyArray_DATA(dirmat_AO);
    printf("state 4\n");

    Fmat_ptr[0]=100;
    int ct;
    ct = 0;
    int i, j;
    for (i=0; i<Nr; i++) {
        for (j=0; j<Nc; j++) {
            Fmat_ptr[j*Nc+i]=i+j;
        }
    }

    printf("just got past changing the matrics\n");


    Py_DECREF(dirmat_AO);
    Py_DECREF(Fmat_AO);

    /* Build the output tuple */
    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *nwops_testNumpyDataExchange(PyObject *self, PyObject *args)
{
    int node, away;
    PyObject *lens, *topo;

    /* Parse the input tuple */
    if (!PyArg_ParseTuple(args, "OO", &lens, &topo))
        return NULL;

    /* Interpret the input objects as numpy arrays. */
    PyObject *len_array = PyArray_FROM_OTF(lens, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *topo_array = PyArray_FROM_OTF(topo, NPY_INT, NPY_IN_ARRAY);

    /* If that didn't work, throw an exception. */
    if (len_array == NULL || topo_array == NULL ) {
        Py_XDECREF(len_array);
        Py_XDECREF(topo_array);
        return NULL;
    }

    /* How many data points are there? */
    int n_topo, n_lens;
    n_topo = (int)PyArray_DIM(topo_array, 0);
    n_lens = (int)PyArray_DIM(len_array, 0);
    printf("n_topo: %d,  n_lens: %d", n_topo, n_lens);

    /* Get pointers to the data as C-types. */
    double *len_ptr    = (double*)PyArray_DATA(len_array);
    int *topo_ptr    = (int*)PyArray_DATA(topo_array);

    /* Call the external C function to compute the chi-squared. */
//    double value = maxLengthToTip( len_ptr, topo_ptr, node, away);
    int i;
    double val;
    for (i=0; i<n_topo; i++) {
        topo_ptr[i] = i;
    }
    for (i=0; i<n_lens; i++) {
        len_ptr[i] += 0.5;
    }
    val = len_ptr[n_lens - 1] + 3.14;


    /* Clean up. */
    Py_DECREF(len_array);
    Py_DECREF(topo_array);

    if (val < 0.0) {
        PyErr_SetString(PyExc_RuntimeError,
                        "maxLengthToTip returned an impossible value.");
        return NULL;
    }

    /* Build the output tuple */
    PyObject *ret = Py_BuildValue("d", val);
    return ret;
}