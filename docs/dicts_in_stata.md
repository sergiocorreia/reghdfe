# Dicts in Stata

Dictionaries, also known as associative arrays, maps, or symbol tables have an unusual syntax in Stata, compared to e.g. Python.

## Python Equivalency

| *Python*          | *Stata*                                                              |
|-------------------|----------------------------------------------------------------------|
| A = dict()        | A = asarray_create([keytype,keydim,minsize,minratio,maxratio])       |
| A[key] = a        | asarray(A, key, a)                                                   |
| a = A[key]        | a = asarray(A, key)                                                      |
| del A[key]        | asarray_remove(A, key)                                               |
| key in A          | asarray_contains(A, key)                                             |
| N = len(A)        | asarray_elements(A)                                                  |
| k = A.keys()      | k = asarray_keys(A)                                                  |
| N/A               | loc=_first(A) loc=_next(A,loc) k=*_key(A, loc), a=*_contents(A, loc) |
| a = A.get(k, 999) | asarray_notfound(A, 999)                                             |

## Optimization

- keytype: either string, real, complex, pointer (default is string). Surround by quotes.
- keydim: if keydim=2 and keytype=real, we have a (very slow) sparse matrix
- minsize: initial size of the hash table; default 100.
- minratio: ?
- maxratio: ?

## Iterating over an array

```mata
for (loc=asarray_first(A); loc!=NULL; loc=asarray_next(A, loc)) {
  k = asarray_key(A, loc)
  v = asarray_contents(A, loc)
}
```

## Hash Function

An alternative may be to use `hash1(x[,n])`

## Counter

```mata
// create counter
counter = asarray_create()
asarray_notfound(counter, 0)

// increase counter
name = "John"
asarray(counter, name, asarray(counter, name) + 1)
```
