#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdint.h>
#include <math.h>
#include <cairo-svg.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>


typedef struct{
  char * lFileName;
  char * rFileName;
  char * oFileName;
  int side;
  int nlayers;
} popt;


void popt_init(popt * p)
{
  // Initialize/set to default
  p->lFileName = NULL;
  p->rFileName = NULL;
  p->oFileName = malloc(1024*sizeof(char));
  sprintf(p->oFileName, "pizza.svg");
  p->side = 101;
  p->nlayers = 5;
  return;
}

/* TODO:
 *  - Don't initiate chrs with 0 pixels or allow them to disappear.
 *  * Determine what order the chrs should be placed in order to avoid twists
 *   (hint, look at cumsum balance of the table).
 *
 */

uint8_t cmap[] = {255,255,255, 	
  240,163,255,
  0,117,220, 	
  153,63,0, 	
  76,0,92, 	
  25,25,25, 	
  0,92,49, 	
  43,206,72, 	
  255,204,153, 	
  128,128,128, 	
  148,255,181, 	
  143,124,0, 	
  157,204,0, 	
  194,0,136, 	
  0,51,128, 	
  255,164,5, 	
  255,168,187, 	
  66,102,0, 	
  255,0,16, 	
  94,241,242, 	
  0,153,143, 	
  224,255,102, 	
  116,10,255, 	
  153,0,0, 	
  255,255,128, 	
  255,255,0, 	
  255,80,5};


void create_svg(double * M, size_t m, size_t n, char * fname)
{
  double f = 10; // scaling factor
  cairo_surface_t * surface = (cairo_surface_t*) cairo_svg_surface_create(fname, f*m, f*n);
  assert(surface != NULL);

  cairo_t * cr = cairo_create(surface);
  cairo_set_antialias (cr, CAIRO_ANTIALIAS_SUBPIXEL); // NOTE: not supported by any backend
  assert(cr != NULL);
  

  double swell = 0.07;
  cairo_set_line_width(cr, 0.0);

  for(int mm = 0; mm<m; mm++)
  {
    for(int nn = 0; nn<n; nn++)
    {
      int L = M[mm + nn*n];
      assert(L>-1);
      assert(L<25);
      double r = cmap[3*L];
      double g = cmap[3*L+1];
      double b = cmap[3*L+2];
      //      cairo_set_line_width(cr, 0);
      cairo_set_source_rgb(cr, r/255, g/255, b/255);
      cairo_rectangle(cr, f*(mm-swell), f*(nn-swell), f*(1.0+swell), f*(1.0+swell));
      cairo_fill(cr);
    }
  }

  cairo_destroy(cr);
  cairo_surface_destroy(surface);

  return;
}

void create_sector_svg(double * T, int * Order, int nrow, int ncol , char * fname)
{
  double width = 501;
  double mid = width/2.0;


  cairo_surface_t * surface = (cairo_surface_t*) cairo_svg_surface_create(fname, width, width);
  assert(surface != NULL);
  cairo_t * cr = cairo_create(surface);
  cairo_set_line_width (cr, 5);
  assert(cr != NULL);

 
  for(int rr = 0; rr<nrow; rr++) // From centre to periphery
  {
  double ang_offset = 0;
  printf("rr = %d\n", rr);
  for(int cc = 0; cc<ncol-1; cc++)
  {
    int L = Order[cc];

    double r = cmap[3*L];
    double g = cmap[3*L+1];
    double b = cmap[3*L+2];

    printf("L = %d, rgb: %f %f %f", L, r, g, b);

    double r0 = mid*rr/nrow;
    double r1 = mid*(rr+1)/nrow;
    double ang0 = ang_offset;
      ang_offset+=(T[L*nrow+rr]*2*M_PI);
    double ang1 = ang_offset;
    printf(" [%f, %f]\n", ang0, ang1);
    assert(ang1>=ang0);
    assert(ang1-2*M_PI < 0.01);
  

    double x0 = r0*sin(ang0)+mid;
    double y0 = r0*cos(ang0)+mid;

    //double x1 = r1*sin(ang1)+mid;
    //double y1 = r1*cos(ang1)+mid;

 //   printf("p0 = (%f, %f), p1 = (%f, %f)\n", x0, y0, x1, y1);


    cairo_set_source_rgb(cr, r/255, g/255, b/255);
    // Draw sector from radius r0 to r1, 
    //             from angle ang0 to ang 1
    cairo_move_to(cr, x0, y0);
    cairo_new_path(cr);
    cairo_arc(cr, mid, mid, r0, ang0, ang1);
    //cairo_line_to(cr, x1, y1);
    cairo_arc_negative(cr, mid, mid, r1, ang1, ang0);
    cairo_close_path(cr);
    //cairo_stroke(cr); 
    cairo_fill(cr);
  }
  }

  cairo_destroy(cr);
  cairo_surface_destroy(surface);
  return;
}


void showTable(FILE * f, double * T, size_t nrow, size_t ncol)
{
  // Show table for first 5 chromosomes
  for(size_t rr = 0; rr<nrow; rr++)
  {
    for(size_t cc = 0; cc<ncol; cc++)
    {
      fprintf(f, "%3.2f\t", T[rr + cc*nrow]);
    }
    fprintf(f, "\n");
  }
  return;
}


uint8_t * readL(char * fname, size_t * NEL)
{
  /* Read label matrix pointed to by p->lfname
  */

  fprintf(stdout, "Reading L-labels from %s\n", fname);

  FILE * f = fopen(fname, "rb");
  fseek(f, 0, SEEK_END); // seek to end of file
  size_t fsize = ftell(f); // get current file pointer
  fseek(f, 0, SEEK_SET); // seek back to beginning of file
  size_t N = fsize/sizeof(uint8_t);

  uint8_t * L = malloc(N*sizeof(uint8_t));


  size_t status = fread(L, sizeof(uint8_t), N, f);    
  fclose(f);

  if(status != N)
  {
    printf("Problems at %d\n", __LINE__);
  }

  NEL[0] = N;
  return L;
}


double * readG(char * fname, size_t * nread)
{
  /* Read label matrix pointed to by p->lfname
  */

  fprintf(stdout, "Reading R from %s\n", fname);

  FILE * f = fopen(fname, "rb");
  fseek(f, 0, SEEK_END); // seek to end of file
  size_t fsize = ftell(f); // get current file pointer
  fseek(f, 0, SEEK_SET); // seek back to beginning of file
  size_t N = fsize/sizeof(double);

  double * G = malloc(N*sizeof(double));

  size_t status = fread(G, sizeof(double), N, f);    
  nread[0] = status;
  fclose(f);

  if(status != N)
  {
    printf("Problems at %d\n", __LINE__);
  }

  return G;
}


void normalizeRows(double * T, size_t nrow, size_t ncol)
{
  // Normalize each row to 1.
  for(size_t rr = 0; rr<nrow; rr++)
  {
    double rsum = 0;
    for(size_t cc = 0; cc<ncol; cc++)
    {
      rsum += T[cc*nrow + rr];
    }

    for(size_t cc = 0; cc<ncol; cc++)
    {
      T[cc*nrow + rr]/=rsum;
    }
  }
}

int radius_to_layer(double r, int nrow)
{
  int layer = round(r*(double) nrow - .5);
  if(layer < 0)
    layer = 0;
  if(layer == nrow)
    layer = nrow-1;
  //   printf("r: %f l: %d/%d\n", r, layer, nrow-1);
  assert(layer>=0);
  assert(layer<nrow);
  return layer;
}

int topo_ok(double * I, int M, size_t pos)
{
  int c = 0;
  int L = (int) I[pos];
  c += (I[pos+1  ] == L) & (I[pos+1+M] != L); // +M
  c += (I[pos+1+M] == L) & (I[pos  +M] != L); // -1
  c += (I[pos  +M] == L) & (I[pos-1+M] != L); // -1
  c += (I[pos-1+M] == L) & (I[pos-1  ] != L); // -M
  c += (I[pos-1  ] == L) & (I[pos-1-M] != L); // -M
  c += (I[pos-1-M] == L) & (I[pos  -M] != L); // +1
  c += (I[pos  -M] == L) & (I[pos+1-M] != L); // +1
  c += (I[pos+1-M] == L) & (I[pos+1  ] != L); // +M

  c += (I[pos+1  ] != L) & (I[pos+1+M] == L); // +M
  c += (I[pos+1+M] != L) & (I[pos  +M] == L); // -1
  c += (I[pos  +M] != L) & (I[pos-1+M] == L); // -1
  c += (I[pos-1+M] != L) & (I[pos-1  ] == L); // -M
  c += (I[pos-1  ] != L) & (I[pos-1-M] == L); // -M
  c += (I[pos-1-M] != L) & (I[pos  -M] == L); // +1
  c += (I[pos  -M] != L) & (I[pos+1-M] == L); // +1
  c += (I[pos+1-M] != L) & (I[pos+1  ] == L); // +M

  int sum = 0;
  for(int aa = -1; aa<2; aa++)
  {
    for(int bb = -1; bb<2; bb++)
    {
      if(I[pos+aa+M*bb] == L)
      {
        sum++;
      }
    }
  }

  if((sum==8)  || (c != 2))
  {
    return 0;
  }

  return 1;
}

int  most_probable(double * T, double * TI, size_t nrow, size_t ncol, 
    int row,
    double p1, 
    double p2,
    double p3,
    double p4)
{

  // 1 Figure out most deficient label for the given radius
  int chr = 0;
  double max_de = 1e99;
  //  printf("row: %d\n", row);
  for(int cc = 1; cc<ncol; cc++)
  {
    double e0 =  pow(TI[cc*nrow + row]-T[cc*nrow + row], 2);
    double e1 =  pow(TI[cc*nrow + row]+1-T[cc*nrow + row], 2);
    double de = e1 - e0;
    //    printf("cc: %d e0: %f - e1: %f  = %f\n", cc, e0, e1, de);
    if(de < max_de)
    {
      max_de = de;
      chr = cc;
    }
  }

  //printf("chr: %d\n", chr);
  //exit(0);
  // 2 If any neighbour has that label, proceed
  // 3 update TI
  return chr;
}

void update_TI(double * TI, double * I, double * RI, int nrow, int ncol, int M, int N)
{

  memset(TI, 0, nrow*ncol*sizeof(double));

  for(size_t kk = 0; kk<M*N; kk++)
  {
    if(RI[kk] <= 1)
    {    
      int col = (int) I[kk];
      int row = radius_to_layer(RI[kk], nrow);
      assert(col<ncol);
      TI[col*nrow + row]++;
      TI[row]++; // count all pixels of the layer
    }
  }
  return;
}
double getTerror(double * T, double * TI, int nrow, int ncol)
{
  double err = 0;
  for(int rr = 0; rr<nrow; rr++)
    for(int cc = 1; cc<ncol; cc++)
    {
      err += pow(T[rr+nrow*cc] - TI[rr+nrow*cc], 2);
    }
  return(err);
}


void eat_tails(double * I, double * RI, int M, double * TI, int nrow, int ncol)
{

  for(size_t pos = 0; pos<M*M; pos++)
  {
    if(RI[pos]<1)
    {
      int label1 = I[pos];
      int sum = 0;
      int label2 = 0;
      for(int aa = -1; aa<2; aa++)
      {
        for(int bb= -1; bb<2; bb++)
        {
          if((abs(aa) + abs(bb)) == 1)
          {
            if(I[pos+aa+bb*M] == label1)
            {
              sum++;
            } else {
              if(label2 == 0)
              {
                label2 = I[pos+aa+bb*M];
              }
            }
          }
        }
      }
      if(label2 > 0)
      {
        if(sum == 1) // End of tail.
        {
          int layer = radius_to_layer(RI[pos], nrow);
          I[pos] = label2;
          //     printf("%zu : %d -> %d\n", pos, label1, label2);
          TI[layer + nrow*label1]--;
          TI[layer + nrow*label2]++;
        }
      }
    }
  }
  return;
}


double optimize_at_borders(double * I, // image
    double * RI, // pre-computed radius
    size_t M, size_t N, // image size
    double * T, // Table with optimal percentage of pixels per layer and class
    size_t nrow, size_t ncol)
{
  // Create a table that holds the current probabilities
  double * TI = malloc(nrow*ncol*sizeof(double));
  update_TI(TI, I, RI, nrow, ncol, M, N);


  printf("TI -- initial:\n");
  showTable(stdout, TI, nrow, ncol);

  // Normalize probability table to number of pixels
  for(int rr = 0; rr<nrow; rr++)
  {
    for(int cc = 1; cc<ncol; cc++)
    {
      T[rr+nrow*cc] *= TI[rr];
    }
  }

  printf("T -- as expected number of pixels:\n");
  showTable(stdout, T, nrow, ncol);
  char  tFileName[] = "pizza.svg.wanted.tsv";
  printf("Writing expected number of pixels per layer and chromosome to '%s'\n", tFileName);
  FILE * tFile = fopen(tFileName, "w");
  showTable(tFile, T, nrow, ncol);
  fclose(tFile);


  double errorOld = -1;
  size_t nIter = 1000;
  for(size_t kk = 0 ; kk < nIter; kk++)
  {
    //for(size_t pos = M+1 ; pos < M*N-M-1; pos++)
    for(size_t pp = M+1 ; pp < 2*M*N; pp++)
    {

      size_t pos = rand() % (M*N);

      if(RI[pos] < 1)
      {
        if(topo_ok(I, M, pos)) // Changing the label at pos does not change the topology
        {

          // Find a direction to pick the new label from
          int step = 0;
          for(int aa = -1; aa<2; aa++)
          {
            for(int bb = -1; bb<2; bb++)
            {
              if(abs(aa)+abs(bb) == 1)
              {
                if((I[pos+aa+M*bb] != I[pos]) & (I[pos+aa+M*bb] != 0))
                {
                  step = aa+M*bb;
                }
              }
            }
          }


          if(step != 0)
          {

            int labelP = I[pos]; // Previous/current label
            int labelN = I[pos+step]; // Label to potentially switch to
            double radius = RI[pos];
            int layer = radius_to_layer(radius, nrow);


            // Try: is energy lowered by switching I[pos] to labelN ?
            double e0 =  pow(TI[labelP*nrow + layer]    -T[labelP*nrow + layer], 2)
              + pow(TI[labelN*nrow + layer]    -T[labelN*nrow + layer], 2);
            double e1 =  pow((TI[labelP*nrow + layer]-1.0)-T[labelP*nrow + layer], 2)
              + pow((TI[labelN*nrow + layer]+1.0)-T[labelN*nrow + layer], 2);

            double de = e1 - e0;

            int use  = 0; 
            if(de<0)
              use = 1;

            if(use == 1){
              double error0 = getTerror(T, TI, nrow, ncol);
              TI[labelN*nrow + (int) layer] ++;
              TI[labelP*nrow + (int) layer] --;
              I[pos] = labelN;
              double error1 = getTerror(T, TI, nrow, ncol);
              assert(error1<=error0);
              assert(fabs( (e1-e0)  - (error1 -error0)) < 0.0001);

            }
          }
        }
      }
    }
    // TODO: Don't do this the last 10 iterations
    //    if(( (int) nIter- (int) kk ) > 10)
    {
      eat_tails(I, RI, M, TI, nrow, ncol);
    }
    double errorNew = getTerror(T, TI, nrow, ncol);
    if( kk % 100 == 0){
      printf("Error at %zu: %f\n", kk,  errorNew);
    }
    if(errorNew == errorOld)
      kk = 999999999999999;
    errorOld = errorNew;

  }

  printf("TI -- final:\n");
  showTable(stdout, TI, nrow, ncol);
  char  tiFileName[] = "pizza.svg.final.tsv";
  printf("Writing final number of pixels per layer and chromosome to '%s'\n", tiFileName);
  FILE * tiFile = fopen(tiFileName, "w");
  showTable(tiFile, TI, nrow, ncol);
  fclose(tiFile);

  if(0){
    printf("TI -- again -- verify that they are the same:\n");
    update_TI(TI, I, RI, nrow, ncol, M, N);
    showTable(stdout, TI, nrow, ncol);
  }

  double err = getTerror(T, TI, nrow, ncol);

  free(TI);
  return err;
}

double optimize(double * I, double * RI, size_t M, size_t N, double * T, size_t nrow, size_t ncol)
{

  // Create a table that holds the current probabilities
  double * TI = malloc(nrow*ncol*sizeof(double));
  memset(TI, 0, nrow*ncol*sizeof(double));
  size_t n2 = 0;
  for(size_t kk = 0 ; kk<M*N; kk++)
  {
    if(I[kk] == 2)
      n2++;
  }
  printf("%zu 2s\n", n2);

  for(size_t kk = 0; kk<M*N; kk++)
  {
    if(RI[kk] <= 1)
    {    
      int col = (int) I[kk];
      int row = radius_to_layer(RI[kk], nrow);
      assert(col<ncol);
      TI[col*nrow + row]++;
    }
  }


  printf("TI -- initial:\n");
  showTable(stdout, TI, nrow, ncol);

  // Normalize probability table to number of pixels
  for(int rr = 0; rr<nrow; rr++)
  {
    for(int cc = 1; cc<ncol; cc++)
    {
      T[rr+nrow*cc] *= TI[rr];
    }
  }

  printf("T -- as expected number of pixels:\n");
  showTable(stdout, T, nrow, ncol);

  for(size_t kk = 0 ; kk < 1000; kk++)
  {
    for(size_t pos = 0 ; pos < M*N; pos++)
    {
      if(I[pos] == 0 && RI[pos] < 1) // if unset
      {      
        double radius = RI[pos];
        double layer = radius_to_layer(radius, nrow);

        int label = most_probable(T, TI, nrow, ncol, 
            layer,
            I[pos + 1], 
            I[pos - 1],
            I[pos + M],
            I[pos - M]);
        int use = 0;
        assert(label>0);

        // printf("r: %f layer: %f, mp: %d\n", radius, layer, label);

        // Check if any of the 4-connected neighbours has this label        
        if(I[pos + 1] == label)
          use = 1;
        if(I[pos -1] == label)
          use = 1;
        if(I[pos + M] == label)
          use = 1;
        if(I[pos - M] == label)
          use = 1;
        use = 1;
        if(use == 1){
          TI[label*nrow + (int) layer] ++;
          I[pos] = label;
        }
      }
    }
  }

  printf("TI -- final:\n");
  showTable(stdout, TI, nrow, ncol);

  double err = 0;
  for(int rr = 0; rr<nrow; rr++)
    for(int cc = 1; cc<ncol; cc++)
    {
      err += pow(T[rr+nrow*cc] - TI[rr+nrow*cc], 2);
    }

  free(TI);
  return err;
}

void initialize_I_random(double * I, double * RI, size_t M, size_t N, int ncol)
{
  int * hasLabel = calloc(ncol, sizeof(int));

  int nSet = 0;
  for(size_t kk = 0; kk<M*N; kk++)
  {
    I[kk] = 0;
    if(RI[kk] < 1)
    {
      if(rand() < 0.01*RAND_MAX)
      {
        int label = (rand()%(ncol)) + 1;
        assert(label>0);
        assert(label<ncol);

        if(hasLabel[label] == 0)
        {
          I[kk] = label;
          nSet++;
          hasLabel[label] = 1;
        }
      }
    }
  }
  printf("%d pixels set\n", nSet);

  return;
}


void initialize_I_pizza(double * I, double * RI, size_t M, size_t N, int ncol, int * order)
{
  int s = (M-1)/2;
  assert(M == N);
  for(int mm = 0; mm< M; mm++)
  {
    for(int nn= 0; nn< M; nn++)
    {
      size_t pos = mm+M*nn;
      if(RI[pos] < 1)
      {
        double x = mm -s;
        double y = nn -s;
        double angle = atan2(x,y);

        int labelIdx = floor((1+angle/M_PI)/2*(ncol-1));
        if(labelIdx>(ncol-2))
          labelIdx = ncol-2;
        int label = order[labelIdx];
        //        printf("labelIdx: %d, label: %d\n", labelIdx, label);
        assert(label>0);
        assert(label<=ncol);
        I[pos] = label;
        if(RI[pos]<0.10)
          I[pos] = 1;
      } else {
        I[pos] = 0;
      }

    }
  }
  return;
}


void topo_ok_test()
{
  double * I = malloc(9*sizeof(double));

  I[0] = 0; I[3] = 0;  I[6] = 0;
  I[1] = 0; I[4] = 1;  I[7] = 1;
  I[2] = 1; I[5] = 1;  I[8] = 1;
  assert(topo_ok(I, 3, 4) == 1);

  I[0] = 0; I[3] = 0;  I[6] = 0;
  I[1] = 0; I[4] = 1;  I[7] = 1;
  I[2] = 1; I[5] = 0;  I[8] = 1;
  assert(topo_ok(I, 3, 4) == 0);

  I[0] = 1; I[3] = 0;  I[6] = 1;
  I[1] = 1; I[4] = 1;  I[7] = 1;
  I[2] = 1; I[5] = 0;  I[8] = 1;
  assert(topo_ok(I, 3, 4) == 0);

  I[0] = 1; I[3] = 0;  I[6] = 1;
  I[1] = 1; I[4] = 1;  I[7] = 1;
  I[2] = 1; I[5] = 1;  I[8] = 1;
  assert(topo_ok(I, 3, 4) == 0);

  // Tricky case: At edge this is fine, but it can possibly create "lines"
  I[0] = 0; I[3] = 0;  I[6] = 0;
  I[1] = 0; I[4] = 1;  I[7] = 1;
  I[2] = 0; I[5] = 1;  I[8] = 1;
  assert(topo_ok(I, 3, 4) == 1);

  I[0] = 0; I[3] = 0;  I[6] = 0;
  I[1] = 0; I[4] = 1;  I[7] = 1;
  I[2] = 0; I[5] = 1;  I[8] = 0;
  assert(topo_ok(I, 3, 4) == 0);

  I[0] = 0; I[3] = 0;  I[6] = 0;
  I[1] = 0; I[4] = 1;  I[7] = 1;
  I[2] = 1; I[5] = 1;  I[8] = 0;
  assert(topo_ok(I, 3, 4) == 0);



  free(I);

  return;
}

double order_error(double * T, int nrow, int ncol)
{
  double * TP = malloc(nrow*ncol*sizeof(double));
  memcpy(TP, T, nrow*ncol*sizeof(double));

  int v= 0;
  for(int rr = 0; rr<nrow; rr++)
  {
    TP[rr] = 0;
    for(int cc = 1; cc<ncol; cc++)
    {
      TP[rr+cc*nrow] += TP[rr+(cc-1)*nrow];
    }
  }
  if(v)
  {
    printf("Column errors:\n");
  }
  double error = 0;
  for(int cc = 1; cc<ncol; cc++)
  {
    double min = 1; double max = 1;
    for(int rr = 0; rr<nrow; rr++)
    {
      double v = TP[rr+cc*nrow];
      if(v<min)
        min = v;
      if(v>max)
        max = v;
    }
    double cerror = pow(max-min,2);

    error += cerror;
    if(v) {
      printf("%.2f ", cerror); }
  }
  if(v) {
    printf(" = %f\n", error); }
  free(TP); 
  return error;
}

void order_reorder(double * TP, double * T, int nrow, int ncol, int * order)
{
  memcpy(TP, T, nrow*sizeof(double));
  for(int cc = 1; cc<ncol; cc++)
  {
    int ccfrom = order[cc-1];
    memcpy(TP+cc*nrow, T+ccfrom*nrow, nrow*sizeof(double));
  }
  return;
}
void int_swap(int * D, size_t n)
{ // swap two element of vector D of size n
  size_t a = rand() % n;
  size_t b = a;
  while(b == a)
  {
    b = rand() % n;
  }
  int t = D[a];
  D[a] = D[b];
  D[b] = t;
  return;
}

int * suggestedOrder(double * T, int nrow, int ncol)
{

  printf("suggestOrder -- T:\n");
  showTable(stdout, T, nrow, ncol);
  double * TP = malloc(nrow*ncol*sizeof(double));

  // TODO: Loop over different permutations of the columns of T
  memcpy(TP, T, nrow*ncol*sizeof(double));

  int * order = malloc((ncol-1)*sizeof(int));
  for(int kk = 0; kk<ncol-1; kk++)
  {
    order[kk] = kk+1;
  }

  int * order0 = malloc((ncol-1)*sizeof(int));
  memcpy(order0, order, (ncol-1)*sizeof(int));

  printf("Pre-err\n");
  showTable(stdout, TP, nrow, ncol);
  double err0 = order_error(TP, nrow, ncol);

  printf("Initial order: %f\n", err0);
  for(size_t tt = 0; tt<20000 ; tt++)
  {
    order_reorder(TP, T, nrow, ncol, order);
    double err = order_error(TP, nrow, ncol);
    if(err<err0)
    {
      //      printf("Smallest error: %f\n", err);
      err0 = err;
      memcpy(order0, order, (ncol-1)*sizeof(int));
    } else {
      memcpy(order, order0, (ncol-1)*sizeof(int));
    }
    int_swap(order, ncol-1);
  }

  printf("Table shown with the suggestOrder -- TP:\n");
  showTable(stdout, TP, nrow, ncol);
  free(TP);

  printf("Placement order of seeds: ");
  for(int kk = 0; kk<(ncol-1); kk++)
  {
    printf("%d ", order0[kk]);
  }
  printf("\n");

  // TODO: XXX ??? 
  return order0; // or order?
}

void usage()
{
  printf("Usage:\n");
  printf("  --lFile file, \n    set file with chromosome labels, should be of type 'uint8'\n");
  printf("  --rFile file, \n    set file with radial information, should be of type 'double'\n");
  printf("  --size N, \n    set output image size (output image will be NxN)\n");
  printf("  --layers N, \n    specify how many radial layers to use\n");
  return;
}
int argparsing(popt * p, int argc, char ** argv)
{

  struct option longopts[] = {
    { "help",         no_argument,       NULL,   'h' },
    // Data
    { "rFile",        required_argument, NULL,   'R' },
    { "lFile",        required_argument, NULL,   'L' },
    // Settings
    { "size",         required_argument, NULL,   's' },
    { "layers",     required_argument, NULL,   'l' },
    // Settings / Forces
    { NULL,           0,                 NULL,   0   }
  };


  int ch;
  while((ch = getopt_long(argc, argv, "hL:R:s:l:", longopts, NULL)) != -1)
  {
    switch(ch)
    {
      case 'L':
        p->lFileName = malloc(strlen(optarg)+1);
        strcpy(p->lFileName, optarg);
        break;
      case 'R':
        p->rFileName = malloc(strlen(optarg)+1);
        strcpy(p->rFileName, optarg);
        break;
      case 's':
        p->side = atoi(optarg);
        break;
      case 'l':
        p->nlayers = atoi(optarg);
        break;
      case 'h':
        usage();
        exit(0);
        break;
      default:
        printf("Unknown argument\n");
        exit(1);
    }
  }

  if(p->lFileName == NULL)
  {
    fprintf(stderr, "No --lFile given (Hint: Use a L*.uint8 file, same as chromflock).\n");
    exit(1);
  }

  if(p->rFileName == NULL)
  {
    fprintf(stderr, "No --rFileName given\n");
    exit(1);
  }

  return 0;
}


int main(int argc, char ** argv)
{
  // Settings
  popt * opt = malloc(sizeof(popt));
  popt_init(opt);

  argparsing(opt, argc, argv);
  int M = opt->side;
  int N = opt->side;

  topo_ok_test();

  srand(time(NULL)*getpid());
  size_t nL = 0;
  printf("Reading %s\n", opt->lFileName);
  uint8_t * L = readL(opt->lFileName, &nL);
  printf(" Found %zu values\n", nL);

  size_t nG;
  double * R = readG(opt->rFileName, &nG);
  printf(" Found %zu values\n", nG);

  if(nL != nG)
  {
    fprintf(stderr, "Number of read numbers from %s and %s does not match!\n", opt->rFileName, opt->lFileName);
    exit(1);
  }


  size_t ncol = 24; // Number of chromosomes + 1 
  size_t nrow = opt->nlayers; // Number of radial regions

  size_t nT = ncol*nrow;
  double * T = malloc(nT*sizeof(double));
  for(size_t kk = 0 ; kk<nT; kk++)
  {
    T[kk] = 0;
  }

  // Create table
  for(size_t kk = 0 ; kk<nL; kk++)
  {
    int cc = ((int) L[kk]);
    if(cc< ncol)
    {
      /* // If reading raw GPSeq -- don't do that
      R[kk] = 2-R[kk];
      if(R[kk]< 0) 
        R[kk] = 0;
      if(R[kk]>1)
        R[kk] = 1;
        */
      if(isfinite(R[kk]))
      {
        int rr = radius_to_layer(R[kk], nrow);

        if (rr >= nrow)
        {    
          printf("R[%zu] = %f -> %d\n", kk, R[kk], rr);
          rr = nrow-1;
        }
        assert(rr<nrow);
        if(isfinite(R[kk]))
        {
          T[cc*nrow + rr]++;
        }
      }
    }
  }
  printf("Occupancy -- row: layer, column: chromosome\n");
  showTable(stdout, T, nrow, ncol);

  normalizeRows(T, nrow, ncol);

  int * order = suggestedOrder(T, nrow, ncol);

  printf("Table after normalization:\n");
  showTable(stdout, T, nrow, ncol);

  create_sector_svg(T, order, nrow, ncol, "sectors.svg");

  double * I = malloc(M*N*sizeof(double));
  double * RI = malloc(M*N*sizeof(double));

  for(size_t mm = 0; mm<M; mm++)
  {
    for(size_t nn = 0; nn<N; nn++)
    {
      RI[mm + M*nn] = sqrt( pow((double) mm-(double) (M-1)/2, 2) + pow((double) nn -(double) (N-1)/2, 2))  / ((double) ((M-1)/2-1));
      //      printf("R[%zu] = %f\n", mm + M*nn, RI[mm + M*nn]);
    }
  }

  printf("Initializing the [%dx%d] image with slices\n", opt->side, opt->side);
  // initialize_I_random(I, RI, M, N, ncol);
  initialize_I_pizza(I, RI, M, N, ncol, order);


  double err = -1;
  //if(0){
  //  err = optimize(I, RI, M, N, T, nrow, ncol);
  //} else {
    err = optimize_at_borders(I, RI, M, N, T, nrow, ncol);
  //}

  printf("----------------------\n");
  printf("Final error: %f\n", err);

  printf("Writing image to '%s'\n", opt->oFileName);
  create_svg(I, M, N, opt->oFileName);

  if(0){ 
    // Show which pixels that can be switched without breaking topology
    double * TO = malloc(M*N*sizeof(double));
    for(size_t kk = 0; kk<M*N; kk++)
    {
      if(RI[kk]>1)
      {
        TO[kk] = 0;
      } else {
        TO[kk] = topo_ok(I, M, kk)+1;
      }
    }
    create_svg(TO, M, N, "topo.svg");
    free(TO);
  }

  free(RI);
  free(I);
  free(T);
  free(R);
  free(L);
  return 0;
}
