/* allocation and deallocation routines */
/*
  Copyright (C) 213 University of Alberta
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>
/*^*/

#include "myfree.h"

void free1int(int *p)
/*< free1int >*/
{
	free(p);
}

void free2int(int **p)
/*< free2int >*/
{
	free(*p);
	free(p);
}

void free1float(float *p)
/*< free1float >*/
{
	free(p);
}

void free2float(float **p)
/*< free2float >*/
{
	free(*p);
	free(p);
}

void free1complex(sf_complex *p)
/*< free1complex >*/
{
	free(p);
}

void free2complex(sf_complex **p)
/*< free2complex >*/
{
	free(*p);
	free(p);
}

