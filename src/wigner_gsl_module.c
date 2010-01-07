/* -*- coding: utf-8; fill-column: 120 -*-

Copyright (C) 2009 Jochen Küpper <python@jochen-kuepper.de>

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free
Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
*/

#include <Python.h>
#include <gsl/gsl_errno.h>

#include "wigner_gsl.h"

static PyObject *_wigner_drot(PyObject *self, PyObject *args)
{
    int j, m1, m2, status;
    double theta;
    gsl_sf_result result;
    if (!PyArg_ParseTuple(args, "iiid", &j, &m1, &m2, &theta))
        return NULL;
    gsl_set_error_handler_off();
    if(GSL_SUCCESS != (status = gsl_sf_wigner_drot_e(j, m1, m2, theta, &result))) {
        /* error */
        PyErr_SetString(PyExc_RuntimeError, "GSL error in gsl_sf_wigner_drot_e");
    }
    return Py_BuildValue("d", result.val);
}


static PyObject *_wigner_drot_e(PyObject *self, PyObject *args)
{
    int j, m1, m2, status;
    double theta;
    gsl_sf_result result;
    if (!PyArg_ParseTuple(args, "iiid", &j, &m1, &m2, &theta))
        return NULL;
    gsl_set_error_handler_off();
    if(GSL_SUCCESS != (status = gsl_sf_wigner_drot_e(j, m1, m2, theta, &result))) {
        /* error */
        PyErr_SetString(PyExc_RuntimeError, "GSL error in gsl_sf_wigner_drot_e");
    }
    return Py_BuildValue("dd", result.val, result.err);
}


static PyMethodDef _wigner_methods[] = {
    {"drot", _wigner_drot, METH_VARARGS, "Reduced Wigner d rotation matrix."},
    {"drot_e", _wigner_drot_e, METH_VARARGS, "Reduced Wigner d rotation matrix with error bounds."},
    {NULL, NULL, 0, NULL}
};


PyMODINIT_FUNC init_wigner_gsl(void)
{
    (void) Py_InitModule("_wigner_gsl", _wigner_methods);
}
