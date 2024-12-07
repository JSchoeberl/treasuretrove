/*
  sources;

  https://tpenguinltg.wordpress.com/2014/02/15/representing-the-heart-shape-precisely

  https://www.researchgate.net/figure/The-3D-shape-of-a-heart-with-equation-2_fig1_330496498
  
*/  


#include <bla.hpp>
using namespace ngbla;

template <typename T>
T Func (T x, T y, T z)
{
  T t1 = (x*x+ 9.0/4 * y*y + z*z -1);
  return t1*t1*t1 - x*x*z*z*z - 9.0/200 * y*y * z*z*z;
}



Vec<3> DFunc (double x, double y, double z)
{
  AutoDiff<3> adx(x,0);
  AutoDiff<3> ady(y,1);
  AutoDiff<3> adz(z,2);

  auto res = Func(adx, ady, adz);
  return { res.DValue(0), res.DValue(1), res.DValue(2) };
}



void Project (Vec<3> & p)
{
  double f = Func(p(0), p(1), p(2));
  Vec<3> df = DFunc(p(0), p(1), p(2));
  p -= f/InnerProduct(df,df) * df;
}


template <typename FUNC>
void FindLevelSet (Vec<3> points[], double tetvals[], FUNC lambda)
{
  bool same01 = (tetvals[0]>0) == (tetvals[1]>0);
  bool same02 = (tetvals[0]>0) == (tetvals[2]>0);
  bool same03 = (tetvals[0]>0) == (tetvals[3]>0);

  
  if (same01 && same02 && same03) return;

  
  Vec<3> ep[4][4];
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      {
        /*
        // tetvals[i] + lam * (tetvals[j]-testvals[i]) = 0
        double lam = -tetvals[i] / (tetvals[j]-tetvals[i]);
        ep[i][j] = points[i] + lam * (points[j]-points[i]);
        */

        // bisection 
        if ((tetvals[i] > 0) == (tetvals[j]>0))
          ep[i][j] = 0.5*(points[i]+points[j]);
        else
          {
            // p = points[i] + t * (points[j]-points[i]);
            double a = 0, b = 1;
            double fa = tetvals[i];
            // double fb = tetvals[j];
            Vec<3> p;
            while (b-a > 1e-12)
              {
                double c = (a+b)/2;
                p = points[i] + c * (points[j]-points[i]);
                double fc = Func(p(0), p(1), p(2));
                if ( (fa>0) == (fc>0) )
                  {
                    a = c;
                    fa = fc;
                  }
                else
                  {
                    b = c;
                    // fb = fc;
                  }
              }
            ep[i][j] = p;
          }
      }

  if (!same01 && same02 && same03)
    {
      if (tetvals[1] > 0)
        lambda (ep[2][1], ep[3][1], ep[0][1]);
      else
        lambda (ep[3][1], ep[2][1], ep[0][1]);
      return;
    }

  if (same01 && !same02 && same03)
    {
      if (tetvals[2] < 0)      
        lambda (ep[3][2], ep[0][2], ep[1][2]);
      else
        lambda (ep[0][2], ep[3][2], ep[1][2]);
      return;
    }
  if (same01 && same02 && !same03)
    {
      if (tetvals[3] > 0)            
        lambda (ep[0][3], ep[1][3], ep[2][3]);
      else
        lambda (ep[1][3], ep[0][3], ep[2][3]);
      return;
    }
    
  if (!same01 && !same02 && !same03)
    {
      if (tetvals[0] < 0)                  
        lambda (ep[0][1], ep[0][2], ep[0][3]);
      else
        lambda (ep[0][2], ep[0][1], ep[0][3]);
      return;
    }

  if (same01 && !same02 && !same03)
    {
      if (tetvals[0] > 0)
        {
          lambda (ep[0][2], ep[1][2], ep[1][3]);      
          lambda (ep[1][3], ep[0][3], ep[0][2]);
        }
      else
        {
          lambda (ep[1][2], ep[0][2], ep[1][3]);      
          lambda (ep[0][3], ep[1][3], ep[0][2]);
        }
    }

  if (!same01 && same02 && !same03)
    {
      if (tetvals[0] < 0)
        {
          lambda (ep[0][1], ep[2][1], ep[2][3]);      
          lambda (ep[2][3], ep[0][3], ep[0][1]);
        }
      else
        {
          lambda (ep[2][1], ep[0][1], ep[2][3]);      
          lambda (ep[0][3], ep[2][3], ep[0][1]);
        }
    }


  if (!same01 && !same02 && same03)
    {
      if (tetvals[0] > 0)
        {
          lambda (ep[0][1], ep[3][1], ep[3][2]);      
          lambda (ep[3][2], ep[0][2], ep[0][1]);
        }
      else
        {
          lambda (ep[3][1], ep[0][1], ep[3][2]);      
          lambda (ep[0][2], ep[3][2], ep[0][1]);
        }
    }

  
  
}




int main()
{
  int nx = 100, ny = 100, nz=100;
  Tensor<3> vals(nx+1, ny+1, nz+1);
  Tensor<3,Vec<3>> pts(nx+1, ny+1, nz+1);
  
  double xmin=-1.5, xmax=1.5;
  double ymin=-1, ymax=1;
  double zmin=-1.5, zmax=1.5;

  double delta = (ymax-ymin)/ny;
  
  for (int ix = 0; ix <= nx; ix++)
    for (int iy = 0; iy <= ny; iy++)
      for (int iz = 0; iz <= nz; iz++)
        {
          double x = xmin + ix*(xmax-xmin)/nx;
          double y = ymin + iy*(ymax-ymin)/ny;
          double z = zmin + iz*(zmax-zmin)/nz;
          Vec<3> p{x,y,z};
          pts(ix,iy,iz) = p;
          vals(ix,iy,iz) = Func(x,y,z);
          Vec<3> grad = DFunc(x,y,z);

          // if we are too close to the level-set, shift the point away 
          if (fabs (vals(ix,iy,iz)) < 0.01 * L2Norm(grad)*delta)
            {
              if (vals(ix,iy,iz) < 0)
                grad *= -1;
              p += 0.1*delta / L2Norm(grad)*grad;
                
              pts(ix,iy,iz) = p;
              vals(ix,iy,iz) = Func(p(0), p(1), p(2));
            }
          
        }

  /*
    0 ... (0,0,0)
    1 ... (1,0,0)
    2 ... (0,1,0)
    3 ... (1,1,0)
    4 ... (0,0,1)
    5 ... (1,0,1)
    6 ... (0,1,1)
    7 ... (1,1,1)
   */
  static int tets[6][4] =
    {
      { 0, 7, 1, 3 },
      { 0, 7, 3, 2 },
      { 0, 7, 2, 6 },
      { 0, 7, 6, 4 },
      { 0, 7, 4, 5 },
      { 0, 7, 5, 1 }
    };


  ofstream out("heart.stl");
  out.precision(16);
  out << "solid\n";
  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++)
      for (int iz = 0; iz < nz; iz++)
        for (int tetnr = 0; tetnr < 6; tetnr++)
          {
            IVec<3> corners[4];
            for (int j = 0; j < 4; j++)
              {
                corners[j][0] = ix + ((tets[tetnr][j]&1) ? 1 : 0);
                corners[j][1] = iy + ((tets[tetnr][j]&2) ? 1 : 0);
                corners[j][2] = iz + ((tets[tetnr][j]&4) ? 1 : 0);
              }
            
            Vec<3> points[4];
            double tetvals[4];
            for (int i = 0; i < 4; i++)
              {
                points[i] = pts(corners[i][0], corners[i][1], corners[i][2]);
                tetvals[i] = vals(corners[i][0], corners[i][1], corners[i][2]);
              }

            FindLevelSet (points, tetvals, [&] (Vec<3> p1, Vec<3> p2, Vec<3> p3)
            {
              // project points to surface
              Project (p1);
              Project (p2);
              Project (p3);
              
              
              Vec<3> normal = Cross (p2-p1, p3-p1);
              /*
              Vec<3> c = 1.0/3 * (p1+p2+p3);
              Vec<3> grad = DFunc(c(0), c(1), c(2));
              Vec<3> normal = 1/L2Norm(grad)*grad;
              */
              out << "facet normal " << normal(0) << " " << normal(1) << " " << normal(2) << "\n";
              out << "outer loop\n";
              out << "vertex " << p1(0) << " " << p1(1) << " " << p1(2) << "\n";
              out << "vertex " << p2(0) << " " << p2(1) << " " << p2(2) << "\n";
              out << "vertex " << p3(0) << " " << p3(1) << " " << p3(2) << "\n";
              out << "endloop\n";
              out << "endfacet\n";
              
              // write stl triangle
            });
          }

  out << "endsolid\n";
}
