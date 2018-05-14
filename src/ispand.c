/*
   Original date of RCS archived version is  Dec 27  1994
 
 
      $Id: ispand.c,v 1.2 2011/11/18 16:48:50 agnew Exp agnew $
 
      $Log: ispand.c,v $
      Revision 1.2  2011/11/18 16:48:50  agnew
      fixed comment lines

      Revision 1.1  2011/09/23 22:39:01  agnew
      Initial revision
 
 
 Returns the bitwise and of integers n1 and n2
*/
ispand (n1,n2)
int *n1,*n2;
{
  int n3,n4;
  n3 = *n1;
  n4 = *n2;
  return(n3&n4);
}

ispand_(n1,n2)
int *n1,*n2;
{
  int n3,n4;
  n3 = *n1;
  n4 = *n2;
  return(n3&n4);
}
