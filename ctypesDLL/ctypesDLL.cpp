// ctypesDLL.cpp : Defines the exported functions for the DLL application.
//

#include "stdafx.h"
#include <stdio.h>	


// ------------------------------------------------------------------------
// A hello world of C library for Python
// ------------------------------------------------------------------------

#include <Python.h>

//#define EXPORT __declspec(dllexport)
#define EXPORT extern "C" __declspec(dllexport)


EXPORT void hello_world() {
	printf("Hello World\n");
}

EXPORT int add(int a, int b) {
	return a + b;
}

EXPORT int string_len(char *s) {
	return strlen(s);
}

EXPORT const char *is_true(int b) {
	return b ? "TRUE" : "FALSE";
}

EXPORT double div_double(double a, double b) {
	return a / b;
}

EXPORT int inc_by_ref(int a, int *b) {
	*b = a + 1;
	return *b;
}

typedef struct StringObject {
	const char *str;
	int index;
	char c;
	int length;
} StringObject;

EXPORT int check_StringObject(StringObject *o) {
	o->length = strlen(o->str);
	o->c = o->str[o->index];
	return 1;
}

EXPORT int inc_array(long a[], int len) {
	int i;
	for (i = 0; i < len; i++) {
		a[i]++;
	}
	return len;
}



// ------------------------------------------------------------------------
// Performance bake off
// ------------------------------------------------------------------------

EXPORT int first(char *s) {
	return (int)(s[0]);
}


EXPORT void charAt(char *s, int size, int pos, char *ch) {
	*ch = s[pos];
}


typedef struct {
	char *s;
	int  size;
	int  pos;
	char *ch;
} charAt_param;


EXPORT void charAt_struct(charAt_param *param) {
	*(param->ch) = param->s[param->pos];
}


static PyObject *
py_first(PyObject *self, PyObject *args)
{
	const char *str;
	int buf_len;
	if (!PyArg_ParseTuple(args, "s#", &str, &buf_len))
		return NULL;

	int n = (long)str[0];

	PyObject *result = Py_BuildValue("i", n);
	return result;
}


static PyObject *
py_charAt(PyObject *self, PyObject *args)
{
	const char *str;
	int size;
	int pos;
	if (!PyArg_ParseTuple(args, "sii", &str, &size, &pos))
		return NULL;

	char ch = str[pos];

	PyObject *result = Py_BuildValue("c", ch);
	return result;
}



// ------------------------------------------------------------------------
// 2009-07-13:
// Need a dummy Py_InitModule in order to compile this as an Python extension.
// Need to find out how to use distutils to build a Windows shared library that is not a Pythion extension?
//
// Also added some methods for performance comparison
// ------------------------------------------------------------------------

PyMethodDef methods[] = {
	{ "py_first", py_first, METH_VARARGS, ""
	},
	{ "py_charAt", py_charAt, METH_VARARGS, ""
	},
	{ NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC
initmylib()
{
	(void)Py_InitModule("mylib", methods);
}

