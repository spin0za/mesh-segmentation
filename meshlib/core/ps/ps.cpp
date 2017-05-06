#include "ps.h"

using namespace MeshLib::RicciFlow;

//int line_width = 1;
float initial_line_width = 0.3;
#define REAL double
//int  bw_ps = 0;
//int  g_id = 1;

#define XMAX 600
#define XMIN 0
#define YMAX 600
#define YMIN 0

#define XSCALE  280
#define YSCALE  280
#define XCENTER 298
#define YCENTER 298


#define PI 3.14159265358979323846


int count = 0; 

int print_head(const char * fname, FILE ** file, int eps)
{
  printf("Writing %s\n", fname);

  *file = fopen(fname, "w");
 
  if (*file == (FILE *) NULL) {
    printf("  Error:  Could not open %s\n", fname);
    return 1;
  }
  if (eps) {
    fprintf(*file, "%%!PS-Adobe-3.0 EPSF-3.0\n");
  } else {
    fprintf(*file, "%%!PS-Adobe-2.0\n");
  }
  fprintf(*file, "%%%%BoundingBox: %d %d %d %d\n", XMIN,YMIN,XMAX,YMAX);
  fprintf(*file, "%%%%Creator: CCGL\n");
  fprintf(*file, "%%%%EndComments\n\n");
  fprintf(*file, "1 setlinecap\n");
  fprintf(*file, "1 setlinejoin\n");
  fprintf(*file, "%f setlinewidth\n", initial_line_width);
  fprintf(*file, "%d %d moveto\n", XMIN, YMIN);
  fprintf(*file, "%d %d lineto\n", XMAX, XMIN);
  fprintf(*file, "%d %d lineto\n", XMAX, YMAX);
  fprintf(*file, "%d %d lineto\n", XMIN, YMAX);
  fprintf(*file, "closepath\nclip\nnewpath\n");
  return 0;
}


void print_edge_eps( FILE * edgefile, MeshLib::CPoint2  point1, MeshLib::CPoint2  point2, double   xscale,double   yscale, double   xoffset, double   yoffset)
{
      fprintf(edgefile, "%d %d moveto\n",
              (int) ( point1[0] * xscale + xoffset),
              (int) ( point1[1] * yscale + yoffset));
      fprintf(edgefile, "%d %d lineto\nstroke\n",
              (int) ( point2[0] * xscale + xoffset),
              (int) ( point2[1] * yscale + yoffset));
}

//miao
void print_circle_eps(FILE *   circlefile, 
					  MeshLib::CPoint2  center,
					  double   radius,
					  double   xscale, 
				      double   yscale, 
				      double   xoffset, 
				      double   yoffset,
					  double   start,
					  double   end )
{
	fprintf(circlefile,"newpath\n");


		fprintf(circlefile, "%f %f %f %f %f arc\nstroke\n", 
             ( center[0] * xscale + xoffset),
             ( center[1] * yscale + yoffset),
 			 ( radius * xscale), start, end );
	
}

int compute_geodesic_center( MeshLib::CPoint2 start, MeshLib::CPoint2 end , MeshLib::CPoint2 & center)
{
	double m[2][2];
	double v[2];
	double im[2][2];


	m[0][0] = 2 * start[0]; 	m[0][1] = 2 * start[1];
	m[1][0] = 2 * end[0]; 		m[1][1] = 2 * end[1];
	
	v[0] = start[0] * start[0] + start[1] * start[1] + 1.0;
	v[1] = end[0] * end[0] + end[1] * end[1] + 1.0;

	double det = m[0][0] * m[1][1] - m[1][0] * m[0][1];

	if( abs(det) < 1e-5 ) return 0; //straight line

	im[0][0] =  m[1][1]/det;
	im[1][1] =  m[0][0]/det;
	im[0][1] = -m[0][1]/det;
	im[1][0] = -m[1][0]/det;


	center[0] = im[0][0] * v[0] + im[0][1] * v[1];
	center[1] = im[1][0] * v[0] + im[1][1] * v[1];

	return 1;
}

void print_circle( FILE * psfile, MeshLib::CPoint2 start, MeshLib::CPoint2 end )
{
	
	MeshLib::CPoint2 center;

	if(  compute_geodesic_center( start, end,center ) )
	{
		double sa = atan2(start[1]-center[1], start[0]-center[0]);
		double ea = atan2(end[1]-center[1], end[0]-center[0]);

		MeshLib::CPoint2 sv = start - center;
		MeshLib::CPoint2 ev = end   - center;

		double radius = mag( sv );

		//CPoint2 normal = cross( sv, ev );

		double crs = sv[0] * ev[1] - sv[1] * ev[0];

		if( crs < 0 )
		{
			double temp;
			temp = sa;
			sa = ea;
			ea = temp;
		}

		if( sa < 0 )
		{
			sa += 2 * PI;
		}

		if( ea < 0 )
		{
			ea += 2 * PI;
		}

		//miao
		double sd = ( sa *  180.0 / PI );
		double ed = ( ea *  180.0 / PI );

        fprintf(psfile, "%f %f moveto\n",
                 (start[0] * XSCALE + XCENTER),
                 (start[1] * YSCALE + YCENTER));

		print_circle_eps( psfile, center, radius, XSCALE,YSCALE,XCENTER,YCENTER, sd,ed );
 }
	else
	{
		count++;
		print_edge_eps( psfile, start, end , XSCALE,YSCALE, XCENTER,YCENTER );
	}

}




void CPs::print(  const char * psfilename )
{
  FILE *psfile;
  
  if (print_head(psfilename, &psfile,  m_eps)) 
  {
    return;
  }

  //print ideal circle

for( CRFMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  CRicciFlowVertex * pV = *viter;
	  CPoint2 c = pV->huv();
	  c = c - CPoint2(0.5,0.5);
	  c = c * 2.0;
	  pV->huv() = c;
}

  for( CRFMesh::MeshEdgeIterator eiter( m_pMesh); !eiter.end(); eiter ++ )
  {
	  CRicciFlowEdge * edge = *eiter;

	  CRicciFlowVertex * v1 = m_pMesh->edgeVertex1( edge );
	  CRicciFlowVertex * v2 = m_pMesh->edgeVertex2( edge );
		
	  MeshLib::CPoint2  point1 = v1->huv();
	  MeshLib::CPoint2  point2 = v2->huv();
	  //print_edge_eps( psfile, point1, point2, XSCALE, YSCALE, XCENTER, YCENTER );
      //print_circle( psfile, point1, point2 );
  }


  for( CRFMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  CRicciFlowVertex * pV = *viter;

	  MeshLib::CPoint2 c = pV->huv();
	  double  r = exp( pV->u() ) * 2.0;
	  print_circle_eps( psfile, c, r, XSCALE, YSCALE, XCENTER, YCENTER, 0.0, 360.0 );

/*
	  CPoint2 C;
	  double  R;
	  _hyperbolic_euclidean( c, r, C, R );

	  print_circle_eps( psfile, C, R, XSCALE, YSCALE, XCENTER, YCENTER, 0.0, 360.0 );
*/
  }
	

/*			
        fprintf(psfile, "1 0 0 setrgbcolor\n%d setlinewidth\n", 1.1);
		int id = 0;
		id = 1;
		while(id < 44)
		{
		Vertex v1 = mesh.idvertex(3398+id*3700);
		Vertex v2 = mesh.idvertex(3393+id*3700);		
		Point  point1 = mesh.point(v1);
		Point  point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );

		 v1 = mesh.idvertex(3393+id*3700);
		 v2 = mesh.idvertex(3396+id*3700);		
	    point1 = mesh.point(v1);
		point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );

		 v1 = mesh.idvertex(3396+id*3700);
		 v2 = mesh.idvertex(3395+id*3700);		
	    point1 = mesh.point(v1);
		point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );

		 v1 = mesh.idvertex(3395+id*3700);
		 v2 = mesh.idvertex(3394+id*3700);	
	    point1 = mesh.point(v1);
		point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );

		 v1 = mesh.idvertex(3394+id*3700);
		 v2 = mesh.idvertex(3397+id*3700);	
	    point1 = mesh.point(v1);
		point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );

		 v1 = mesh.idvertex(3397+id*3700);
		 v2 = mesh.idvertex(3400+id*3700);
	    point1 = mesh.point(v1);
		point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );

		 v1 = mesh.idvertex(3400+id*3700);
		 v2 = mesh.idvertex(3399+id*3700);
	    point1 = mesh.point(v1);
		point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );

		 v1 = mesh.idvertex(3399+id*3700);
		 v2 = mesh.idvertex(3398+id*3700);
	    point1 = mesh.point(v1);
		point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );
		id = id+6;

		}
*/
  if (!m_eps) {
    fprintf(psfile, "showpage\n");
  }
 
  fprintf(psfile, "0 0 1 setrgbcolor\n");
  print_circle_eps( psfile, CPoint2(0,0), 1.0, XSCALE, YSCALE, XCENTER,YCENTER, 0, 360 ); 

  printf("count  %d\n", count);
  fclose(psfile);
}


void  CPs::_hyperbolic_euclidean( MeshLib::CPoint2 c, double r, MeshLib::CPoint2 & C, double & R ) //converting a hyperbolic circle (c,r) to a Euclidean circle (C,R)
{
	double norm = c[0]*c[0]+c[1]*c[1];
	double mu = (exp(r)-1.0)/(exp(r)+1.0);
	mu = mu * mu;
	double a = 1 - mu * norm;

	double B = 2.0 * ( mu - 1.0) * c[0]/a;
	double E = 2.0 * ( mu - 1.0) * c[1]/a;
	double D = (norm - mu)/a;

	C[0] = -B/2.0;
	C[1] = -E/2.0;
	R    = sqrt( B*B/4.0+E*E/4.0-D);
}




void CPs::hyperbolic_print( const  char * psfilename )
{
  FILE *psfile;
  
  if (print_head(psfilename, &psfile,  m_eps)) 
  {
    return;
  }


  for( CRFMesh::MeshEdgeIterator eiter( m_pMesh); !eiter.end(); eiter ++ )
  {
	  CRicciFlowEdge * edge = *eiter;

	  CRicciFlowVertex * v1 = m_pMesh->edgeVertex1( edge );
	  CRicciFlowVertex * v2 = m_pMesh->edgeVertex2( edge );
		
	  MeshLib::CPoint2  point1 = v1->huv();
	  MeshLib::CPoint2  point2 = v2->huv();
	
      print_circle( psfile, point1, point2 );
  }


  for( CRFMesh::MeshVertexIterator viter( m_pMesh ); !viter.end(); viter ++ )
  {
	  CRicciFlowVertex * pV = *viter;

	  MeshLib::CPoint2 c = pV->huv();
	  double e = exp( pV->u() );
      double r = log(  (1+e)/(1 - e ) );

	  MeshLib::CPoint2 C;
	  double  R;
	  _hyperbolic_euclidean( c, r, C, R );

	  print_circle_eps( psfile, C, R, XSCALE, YSCALE, XCENTER, YCENTER, 0.0, 360.0 );

  }
	

/*			
        fprintf(psfile, "1 0 0 setrgbcolor\n%d setlinewidth\n", 1.1);
		int id = 0;
		id = 1;
		while(id < 44)
		{
		Vertex v1 = mesh.idvertex(3398+id*3700);
		Vertex v2 = mesh.idvertex(3393+id*3700);		
		Point  point1 = mesh.point(v1);
		Point  point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );

		 v1 = mesh.idvertex(3393+id*3700);
		 v2 = mesh.idvertex(3396+id*3700);		
	    point1 = mesh.point(v1);
		point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );

		 v1 = mesh.idvertex(3396+id*3700);
		 v2 = mesh.idvertex(3395+id*3700);		
	    point1 = mesh.point(v1);
		point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );

		 v1 = mesh.idvertex(3395+id*3700);
		 v2 = mesh.idvertex(3394+id*3700);	
	    point1 = mesh.point(v1);
		point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );

		 v1 = mesh.idvertex(3394+id*3700);
		 v2 = mesh.idvertex(3397+id*3700);	
	    point1 = mesh.point(v1);
		point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );

		 v1 = mesh.idvertex(3397+id*3700);
		 v2 = mesh.idvertex(3400+id*3700);
	    point1 = mesh.point(v1);
		point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );

		 v1 = mesh.idvertex(3400+id*3700);
		 v2 = mesh.idvertex(3399+id*3700);
	    point1 = mesh.point(v1);
		point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );

		 v1 = mesh.idvertex(3399+id*3700);
		 v2 = mesh.idvertex(3398+id*3700);
	    point1 = mesh.point(v1);
		point2 = mesh.point(v2);
		print_circle( psfile, point1,point2 );
		id = id+6;

		}
*/
  if (!m_eps) {
    fprintf(psfile, "showpage\n");
  }

  fprintf(psfile, "1 0 0 setrgbcolor\n");
  print_circle_eps( psfile, CPoint2(0,0), 0.005, XSCALE, YSCALE, XCENTER,YCENTER, 0, 360 ); 

  fprintf(psfile, "0 0 1 setrgbcolor\n");
  print_circle_eps( psfile, CPoint2(0,0), 1.0, XSCALE, YSCALE, XCENTER,YCENTER, 0, 360 ); 

  printf("count  %d\n", count);
  fclose(psfile);
}
