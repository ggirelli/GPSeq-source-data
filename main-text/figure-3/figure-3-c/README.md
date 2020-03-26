# Fig 3C

The plot was generated with a custom program. To build it use the supplied makefile, i.e.,

```
make main
```

To generate the image `pizza.svg`, use:

```
./main --lFile L_noy_new3Trans.uint8 --rFile G_noy_new3Trans.double 
```

Unfortunately we didn't think of storing the random seed at the time when the image was created (the code uses `srand(time(NULL)*getpid());`) and a different image will be created each time that the program is run. 

The input to the program consists of 1 Mb bins with an associate chromosome label and radial position. 
 * `L_noy_new3Trans.uint8` describes which chromosome each bin belong to, encoded as unsigned 8-bit integer, i.e., 1 for Chr1, 2 for Chr2, ..., 23 for ChrX.
 * `G_noy_new3Trans.double` says what radial position each bin has, encoded as 64-bit float, i.e., 0 denotes the centre and 1 denotes the lamin.

