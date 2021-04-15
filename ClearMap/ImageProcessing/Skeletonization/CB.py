# -*- coding: utf-8 -*-
"""
Assymetric 3D Thinning Algorithm based on Isthmusses

Reference:
  Asymmetric Parallel 3D Thinning Scheme and Algorithms Based on Isthmuses
  Courpie, Bertrand 2016
"""

import Topology2D as top2d;
import Topology3D as top3d;


def matchK2_z(cube):
  """Matches critical 2-clique in z direction"""
  if not cube[1,1,0] or not cube[1,1,1]:
    return False;
  w = np.logical_or(cube[:,:,0], cube[:,:,1]);
  index = top2d.planeIndex(w);
  indexb = top2d.planeIndex(np.logical_not(w));
  if top2d.t4barLUT[index] == 1 and top2d.t8LUT[indexb] == 1:
    return False;
  else:
    return True;

def matchK2(cube):
  """Matches critical 2-cliques in all 3 directions"""
  if matchK2_z(cube):
    return True;
  cubey = top3d.rotate(cube, axis = 0, steps = 1);
  if matchK2_z(cubey):
    return True;
  cubex = top3d.rotate(cube, axis = 1, steps = 1);
  if matchK2_z(cubex):
    return True;  
  return False;

# ordering
#[0 0 0] [0 0 0] [0 0 0]
#[3 4 5] [12 13 14] [0 0 0]
#[0 1 2] [ 9 10 11] [0 0 0]

#[3 4 5] [12 13 14] [0 0 0]
#[0 1 2] [ 9 10 11] [0 0 0]

def matchK1_z(cube):
  


def matchK1(cube):
  pass

def matchK0_z(cube):  
  pass




static int32_t TopoTab[256][2] = 
{
  {1,0},  {1,1},  {1,1},  {1,1},  {1,1},  {1,1},  {1,1},  {1,1},  /*  0 -  7 */
  {1,1},  {2,2},  {2,2},  {2,2},  {1,1},  {1,1},  {1,1},  {1,1},  /*  8 - 1f */
  {1,1},  {2,2},  {2,2},  {2,2},  {1,1},  {1,1},  {1,1},  {1,1},  /* 10 - 17 */
  {1,1},  {2,2},  {2,2},  {2,2},  {1,1},  {1,1},  {1,1},  {1,1},  /* 18 - 1f */
  {1,1},  {2,2},  {2,2},  {2,2},  {2,2},  {2,2},  {2,2},  {2,2},  /* 20 - 27 */
  {2,2},  {3,3},  {3,3},  {3,3},  {2,2},  {2,2},  {2,2},  {2,2},  /* 28 - 2f */
  {1,1},  {2,2},  {2,2},  {2,2},  {1,1},  {1,1},  {1,1},  {1,1},  /* 30 - 37 */
  {1,1},  {2,2},  {2,2},  {2,2},  {1,1},  {1,1},  {1,1},  {1,1},  /* 38 - 3f */
  {1,1},  {1,1},  {2,2},  {1,1},  {2,2},  {1,1},  {2,2},  {1,1},  /* 40 - 47 */
  {2,2},  {2,2},  {3,3},  {2,2},  {2,2},  {1,1},  {2,2},  {1,1},  /* 48 - 4f */
  {1,1},  {1,1},  {2,2},  {1,1},  {1,1},  {0,1},  {1,1},  {0,1},  /* 50 - 57 */
  {1,1},  {1,1},  {2,2},  {1,1},  {1,1},  {0,1},  {1,1},  {0,1},  /* 58 - 5f */
  {1,1},  {1,1},  {2,2},  {1,1},  {2,2},  {1,1},  {2,2},  {1,1},  /* 60 - 67 */
  {2,2},  {2,2},  {3,3},  {2,2},  {2,2},  {1,1},  {2,2},  {1,1},  /* 68 - 6f */
  {1,1},  {1,1},  {2,2},  {1,1},  {1,1},  {0,1},  {1,1},  {0,1},  /* 70 - 77 */
  {1,1},  {1,1},  {2,2},  {1,1},  {1,1},  {0,1},  {1,1},  {0,1},  /* 78 - 7f */
  {1,1},  {1,1},  {2,2},  {1,1},  {2,2},  {1,1},  {2,2},  {1,1},  /* 80 - 87 */
  {2,2},  {2,2},  {3,3},  {2,2},  {2,2},  {1,1},  {2,2},  {1,1},  /* 88 - 8f */
  {2,2},  {2,2},  {3,3},  {2,2},  {2,2},  {1,1},  {2,2},  {1,1},  /* 90 - 97 */
  {2,2},  {2,2},  {3,3},  {2,2},  {2,2},  {1,1},  {2,2},  {1,1},  /* 98 - 9f */
  {2,2},  {2,2},  {3,3},  {2,2},  {3,3},  {2,2},  {3,3},  {2,2},  /* a0 - a7 */
  {3,3},  {3,3},  {4,4},  {3,3},  {3,3},  {2,2},  {3,3},  {2,2},  /* a8 - af */
  {2,2},  {2,2},  {3,3},  {2,2},  {2,2},  {1,1},  {2,2},  {1,1},  /* b0 - b7 */
  {2,2},  {2,2},  {3,3},  {2,2},  {2,2},  {1,1},  {2,2},  {1,1},  /* b8 - bf */
  {1,1},  {1,1},  {2,2},  {1,1},  {2,2},  {1,1},  {2,2},  {1,1},  /* c0 - c7 */
  {2,2},  {2,2},  {3,3},  {2,2},  {2,2},  {1,1},  {2,2},  {1,1},  /* c8 - cf */
  {1,1},  {1,1},  {2,2},  {1,1},  {1,1},  {0,1},  {1,1},  {0,1},  /* d0 - d7 */
  {1,1},  {1,1},  {2,2},  {1,1},  {1,1},  {0,1},  {1,1},  {0,1},  /* d8 - df */
  {1,1},  {1,1},  {2,2},  {1,1},  {2,2},  {1,1},  {2,2},  {1,1},  /* e0 - e7 */
  {2,2},  {2,2},  {3,3},  {2,2},  {2,2},  {1,1},  {2,2},  {1,1},  /* e8 - ef */
  {1,1},  {1,1},  {2,2},  {1,1},  {1,1},  {0,1},  {1,1},  {0,1},  /* f0 - f7 */
  {1,1},  {1,1},  {2,2},  {1,1},  {1,1},  {0,1},  {1,1},  {0,1}   /* f8 - ff */
};


/* ==================================== */
int32_t t4(int32_t v) /* pour un objet en 4-connexite - v est le masque binaire du voisinage */
/* ==================================== */
{
  v = ~v & 0xff;
  return TopoTab[v][0];
}

/* ==================================== */
int32_t t4b(int32_t v) /* pour un objet en 4-connexite - v est le masque binaire du voisinage */
/* ==================================== */
{
  return TopoTab[v][0];
}




/* ==================================== */
int32_t asym_match_vois2(uint8_t *v)
/* ==================================== */
/*
               12      11      10       
               13       8       9
               14      15      16

		3	2	1			
		4      26	0
		5	6	7
Teste si les conditions suivantes sont r�unies:
1: v[8] et v[26] doivent �tre dans l'objet et simples et non s�lectionn�s
2: for i = 0 to 7 do w[i] = v[i] || v[i+9] ; w[0...7] doit �tre non 2D-simple
Si le test r�ussit, le point 8 est marqu� SELECTED
*/
{
  uint8_t t;
  if (!IS_SIMPLE(v[8]) || !IS_SIMPLE(v[26]) || IS_SELECTED(v[8]) || IS_SELECTED(v[26])) return 0;
  if (v[0] || v[9]) t = 1; else t = 0;
  if (v[1] || v[10]) t |= 2;
  if (v[2] || v[11]) t |= 4;
  if (v[3] || v[12]) t |= 8;
  if (v[4] || v[13]) t |= 16;
  if (v[5] || v[14]) t |= 32;
  if (v[6] || v[15]) t |= 64;
  if (v[7] || v[16]) t |= 128;
  if ((t4b(t) == 1) && (t8(t) == 1)) return 0; // simple 2D
  SET_SELECTED(v[8]);
  return 1;
} // asym_match_vois2()

/* ==================================== */
int32_t asym_match_vois1(uint8_t *v)
/* ==================================== */
// A A  P1 P2  B B
// A A  P3 P4  B B
// avec pour localisations possibles :
// 12 11   3  2   21 20 
// 13  8   4 26   22 17
// et :
// 11 10    2 1   20 19
//  8  9   26 0   17 18
//
// Teste si les trois conditions suivantes sont r�unies:
// 1: (P1 et P4) ou (P2 et P3)
// 2: tous les points Pi non nuls doivent �tre simples et non marqu�s SELECTED
// 3: A et B sont tous nuls ou [au moins un A non nul et au moins un B non nul]
// Si le test r�ussit, un des points Pi non nuls est marqu� SELECTED
{
  int32_t ret = 0;
  if (!((v[2] && v[4]) || (v[3] && v[26]))) goto next1;
  if ((IS_OBJECT(v[2])  && (!IS_SIMPLE(v[2])  || IS_SELECTED(v[2]))) ||
      (IS_OBJECT(v[3])  && (!IS_SIMPLE(v[3])  || IS_SELECTED(v[3]))) ||
      (IS_OBJECT(v[4])  && (!IS_SIMPLE(v[4])  || IS_SELECTED(v[4]))) ||
      (IS_OBJECT(v[26]) && (!IS_SIMPLE(v[26]) || IS_SELECTED(v[26])))) goto next1;
  if ((v[12] || v[11] || v[13] || v[8] || v[21] || v[20] || v[22] || v[17]) &&
      ((!v[12] && !v[11] && !v[13] && !v[8]) || 
       (!v[21] && !v[20] && !v[22] && !v[17]))) goto next1;
  if (v[2])  SET_SELECTED(v[2]);
  else if (v[3])  SET_SELECTED(v[3]);
  else if (v[4])  SET_SELECTED(v[4]);
  else if (v[26]) SET_SELECTED(v[26]);
  ret = 1;
 next1:
  if (!((v[2] && v[0]) || (v[1] && v[26]))) goto next2;
  if ((IS_OBJECT(v[2])  && (!IS_SIMPLE(v[2])  || IS_SELECTED(v[2]))) ||
      (IS_OBJECT(v[1])  && (!IS_SIMPLE(v[1])  || IS_SELECTED(v[1]))) ||
      (IS_OBJECT(v[0])  && (!IS_SIMPLE(v[0])  || IS_SELECTED(v[0]))) ||
      (IS_OBJECT(v[26]) && (!IS_SIMPLE(v[26]) || IS_SELECTED(v[26])))) goto next2;
  if ((v[10] || v[11] || v[9] || v[8] || v[19] || v[20] || v[18] || v[17]) &&
      ((!v[10] && !v[11] && !v[9] && !v[8]) || 
       (!v[19] && !v[20] && !v[18] && !v[17]))) goto next2;
  if (v[2])  SET_SELECTED(v[2]);
  else if (v[1])  SET_SELECTED(v[1]);
  else if (v[0])  SET_SELECTED(v[0]);
  else if (v[26]) SET_SELECTED(v[26]);
  ret = 1;
 next2:
  return ret;
} // asym_match_vois1()

/* ==================================== */
int32_t asym_match_vois0(uint8_t *v)
/* ==================================== */
/*
               12      11      10
               13       8       9
               14      15      16

		3	2	1
		4      26	0
		5	6	7

               21      20      19
               22      17      18
               23      24      25

Teste si les conditions suivantes sont r�unies:
1: au moins un des ensembles {26,12}, {26,10}, {26,14}, {26,21} est inclus dans l'objet, et
2: les points non nuls du cube 2x2x2 contenant cet ensemble sont tous simples, 
   non marqu�s SELECTED
Si le test r�ussit, le point 26 est marqu� SELECTED
*/
{
  if (!v[26]) return 0;
  if (!IS_SIMPLE(v[26]) || IS_SELECTED(v[26])) return 0;
  if (!(v[12] || v[10] || v[14] || v[21])) return 0;
  if (v[12])
  { /*         12      11
               13       8

		3	2
		4      26 */
     if (!IS_SIMPLE(v[12]) || IS_SELECTED(v[12])) return 0;
     if (v[11] && (!IS_SIMPLE(v[11]) || IS_SELECTED(v[11]))) return 0;
     if (v[13] && (!IS_SIMPLE(v[13]) || IS_SELECTED(v[13]))) return 0;
     if (v[ 8] && (!IS_SIMPLE(v[ 8]) || IS_SELECTED(v[ 8]))) return 0;
     if (v[ 3] && (!IS_SIMPLE(v[ 3]) || IS_SELECTED(v[ 3]))) return 0;
     if (v[ 2] && (!IS_SIMPLE(v[ 2]) || IS_SELECTED(v[ 2]))) return 0;
     if (v[ 4] && (!IS_SIMPLE(v[ 4]) || IS_SELECTED(v[ 4]))) return 0;
  }
  if (v[10])
  { /*
               11      10
               8       9

		2	1
		26	0 */
     if (!IS_SIMPLE(v[10]) || IS_SELECTED(v[10])) return 0;
     if (v[11] && (!IS_SIMPLE(v[11]) || IS_SELECTED(v[11]))) return 0;
     if (v[ 8] && (!IS_SIMPLE(v[ 8]) || IS_SELECTED(v[ 8]))) return 0;
     if (v[ 9] && (!IS_SIMPLE(v[ 9]) || IS_SELECTED(v[ 9]))) return 0;
     if (v[ 1] && (!IS_SIMPLE(v[ 1]) || IS_SELECTED(v[ 1]))) return 0;
     if (v[ 2] && (!IS_SIMPLE(v[ 2]) || IS_SELECTED(v[ 2]))) return 0;
     if (v[ 0] && (!IS_SIMPLE(v[ 0]) || IS_SELECTED(v[ 0]))) return 0;
  }
  if (v[14])
  { /*         13       8
               14      15

		4      26
		5	6 */
     if (!IS_SIMPLE(v[14]) || IS_SELECTED(v[14])) return 0;
     if (v[13] && (!IS_SIMPLE(v[13]) || IS_SELECTED(v[13]))) return 0;
     if (v[15] && (!IS_SIMPLE(v[15]) || IS_SELECTED(v[15]))) return 0;
     if (v[ 8] && (!IS_SIMPLE(v[ 8]) || IS_SELECTED(v[ 8]))) return 0;
     if (v[ 6] && (!IS_SIMPLE(v[ 6]) || IS_SELECTED(v[ 6]))) return 0;
     if (v[ 5] && (!IS_SIMPLE(v[ 5]) || IS_SELECTED(v[ 5]))) return 0;
     if (v[ 4] && (!IS_SIMPLE(v[ 4]) || IS_SELECTED(v[ 4]))) return 0;
  }
  if (v[21])
  {  /*		3	2
		4      26

               21      20
               22      17 */
     if (!IS_SIMPLE(v[21]) || IS_SELECTED(v[21])) return 0;
     if (v[17] && (!IS_SIMPLE(v[17]) || IS_SELECTED(v[17]))) return 0;
     if (v[20] && (!IS_SIMPLE(v[20]) || IS_SELECTED(v[20]))) return 0;
     if (v[22] && (!IS_SIMPLE(v[22]) || IS_SELECTED(v[22]))) return 0;
     if (v[ 3] && (!IS_SIMPLE(v[ 3]) || IS_SELECTED(v[ 3]))) return 0;
     if (v[ 2] && (!IS_SIMPLE(v[ 2]) || IS_SELECTED(v[ 2]))) return 0;
     if (v[ 4] && (!IS_SIMPLE(v[ 4]) || IS_SELECTED(v[ 4]))) return 0;
  }
  SET_SELECTED(v[26]);
  return 1;
} // asym_match_vois0()

/* ==================================== */
int32_t asym_match2(uint8_t *v)
/* ==================================== */
{
  int32_t ret = 0;
  if (asym_match_vois2(v)) ret = 1;
  isometrieXZ_vois(v);
  if (asym_match_vois2(v)) ret = 1;
  isometrieXZ_vois(v);
  isometrieYZ_vois(v);
  if (asym_match_vois2(v)) ret = 1;
  isometrieYZ_vois(v);
  return ret;
} /* asym_match2() */

/* ==================================== */
int32_t asym_match1(uint8_t *v)
/* ==================================== */
{
  int32_t ret = 0;
  if (asym_match_vois1(v)) ret = 1;
  isometrieXZ_vois(v);
  if (asym_match_vois1(v)) ret = 1;
  isometrieXZ_vois(v);
  isometrieYZ_vois(v);
  if (asym_match_vois1(v)) ret = 1;
  isometrieYZ_vois(v);
  return ret;
} /* asym_match1() */

/* ==================================== */
int32_t asym_match0(uint8_t *v)
/* ==================================== */
{
  int32_t ret = 0;
  if (asym_match_vois0(v)) ret = 1;
  return ret;
} /* asym_match0() */

/* ==================================== */
int32_t lskelAMK3(struct xvimage *image, 
		   int32_t n_steps,
		   struct xvimage *inhibit)
/* ==================================== */


