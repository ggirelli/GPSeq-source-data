# Fig 3C

The plot was generated with a custom program. To build it use the supplied makefile, i.e.,

```
make main
```

To generate the image `pizza.svg`, use:

```
./main --lFile L_noy_new3Trans.uint8 --rFile G_noy_new3Trans.double 
```

where the data files are the same as for the 10000 structures at 1Mb.

Unfortunately we didn't think of storing the random seed at the time when the image was created (the code uses 'srand(time(NULL)*getpid());')

