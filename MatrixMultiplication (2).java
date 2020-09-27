

/*
Write a program to implement the recursive Strassen's matrix multiplication algorithm in  C++ or Java. Instantiate a matrix
    as an object. Use this matrix multiplication routine to find Eigen values of a matrix by using the QR factorization based method.
 */

/*
Steps to find Eigen value :-
1) Intilizing matrix through Matrix class and Entering data through User(to simplify we have consider matrix order as power of 2 and sqaure matrix).
2) Suppose we denote Matrix by A i.e A=QR-[Q = orthonormal matrix and R is Upper traingular matrix]
3) finding Q using function QR_algorithm_gram_Schmidt() using Gram and schmit algorithm reference.
4) To find eigen value we have to find matrix R. So R=A*Inverse(Q).
5) for finding inverse(Q) we know property ass Q*(Q^T)=I as Q is orthonomal matrix. Hence Inverse(Q)=Q^T.
6) Multiplied Q^T to A to find R using Strassen's recursive algorithm.
7) Now we have R which is diagonal matrix hence to find eigen value we took diagonal value of matrix R.
*/

import java.util.*;

//Matrix class to initilize matrix.
class Matrix 
{
  double c[][];

  Matrix(int n) {
    c = new double[n][n];
  }

}

class Main
{

  public static void main(String[] args) {
    System.out.println("Enter the size of A(it should be in power of 2)");
    Scanner sc = new Scanner(System.in);

    int n = sc.nextInt(); // enter the order of matrix

    Matrix A = new Matrix(n); // Initilizing matrix A using Matrix class.
    enter_matrix(A.c, n); // enter_Matrix function take value from user.
    System.out.println();
    System.out.println("Matrix A: ");
    System.out.println();
    display_matrix(A.c, n); // display matrix entered by user.
    System.out.println();
    Matrix Q = new Matrix(n); // Q is an othro-normal matrix in equation A=Q*R.
    double R[][]=new double[n][n]; // A=QR Here, R is Upper triangular matrix



    // QR_algorithm_gram_Schmidt() is a function which finds Q matrix in equation
    // A=Q*R using gram-schmit algorithm.
    QR_algorithm_gram_Schmidt(A.c, Q.c, n);



    // code to display matrix Q foundd using QR_algorithm_gram_Schmidt function.
    System.out.println("Matrix Q: ");
    System.out.println();
    display_matrix(Q.c, n); // display matrix entered by user.
    System.out.println();


    // as per step 5 to find inverse of Matrix Q whhich is transpose itself.
    
    double QT[][]; 
    // QT matrix stores transpose of Matrix Q which is inverse of Q using property
    // in step 5.2
    QT = transpose(Q.c, n); // Q*Q^t=I
    System.out.println("Inverse of Matrix Q :");
    System.out.println();
    display_matrix(QT, n);
    System.out.println();
  
    // Code for Multiplication between inverse(Q) and A to find R using strassen
    // matrix.
    Strassen_multiplication(QT, A.c, R, 0, 0, 0, 0, n);
    // R is upper-traignular matrix.

    System.out.println("Matrix R:-");
    System.out.println();
    // displaying Matrix R
    display_matrix(R, n);
    System.out.println();


    // As R is upper traingular matrix. Hence by property of upper triangular matrix
    // its Eigen value is diagonal Element.
    display_eigenvalue(R, n); // display eigen value.


    //to check correctness of algorithm, multiplying Matrix Q and Matrix R displaying Matrix A
    double A_after_finding_QR[][]=new double[n][n];
    System.out.println("Matrix A: ");
    Strassen_multiplication(Q.c, R, A_after_finding_QR,0, 0, 0, 0, n);
    System.out.println();
    display_matrix(A_after_finding_QR, n); // display matrix entered by user.
    System.out.println();
  }

  // this function display Eigen Value from R matrix.
  public static void display_eigenvalue(double c[][], int n)
  {
    System.out.println("Eigen values are:");
    System.out.println();
    for (int i = 0; i < n; i++) 
    {
      System.out.println("Eigen Value " + (i + 1) + " " + c[i][i]);
    }
  }

  // This function is used to take user data for matrix.
  public static void enter_matrix(double a[][], int n)
  {
    Scanner sc = new Scanner(System.in);
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        a[i][j] = sc.nextInt();
      }
    }
  }

  // this function finds transpose of matrix.
  public static double[][] transpose(double q[][], int n)
  {
    double transposed_matrix[][] = new double[n][n];
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        transposed_matrix[i][j] = q[j][i];
      }
    }
    return transposed_matrix;
  }

  // display matrix a.
  public static void display_matrix(double a[][], int n) 
  {
    for (int i = 0; i < n; i++) 
    {
      for (int j = 0; j < n; j++) 
      {
        System.out.print(a[i][j] + " ");
      }
      System.out.println();
    }
  }

  // Strassen Multiplication function for matrix a and Matrix b.
  public static void Strassen_multiplication(double a[][], double b[][], double result[][],int i, int j, int k, int l, int n) 
  {

    if (n <= 2) 
    {
      result[i][l] += a[i][j] * b[k][l] + a[i][j + 1] * b[k + 1][l];
      result[i][l+1] += a[i][j] * b[k][l + 1] + a[i][j + 1] * b[k + 1][l + 1];
      result[i+1][l] += a[i + 1][j] * b[k][l] + a[i + 1][j + 1] * b[k + 1][l];
      result[i+1][l+1] += a[i + 1][j] * b[k][l + 1] + a[i + 1][j + 1] * b[k + 1][l + 1];
      
    } else 
    {
      int mid = n / 2;
      
      Strassen_multiplication(a, b ,result, i, j, k, l, mid);
      Strassen_multiplication(a, b,result, i, j + mid, k + mid, l, mid);
      Strassen_multiplication(a, b,result, i, j, k, l + mid, mid);
      Strassen_multiplication(a, b,result, i, j + mid, k + mid, l + mid, mid);
      Strassen_multiplication(a, b,result, i + mid, j, k, l, mid);
      Strassen_multiplication(a, b,result, i + mid, j + mid, k + mid, l, mid);
      Strassen_multiplication(a, b,result, i + mid, j, k, l + mid, mid);
      Strassen_multiplication(a, b,result, i + mid, j + mid, k + mid, l + mid, mid);

    }
  }

  public static void QR_algorithm_gram_Schmidt(double a[][], double q[][], int n) 
  {
    double s[] = new double[n];
    double s1[] = new double[n];

    // s1 stores current column in Q matrix.
    // s array takes all the column before i and helps to stores vector that are
    // perpendicular to i in result array
    // so that they can be remove from vectore current vector under consideration
    // stored in s1 in the end.
    for (int i = 0; i < n; i++) 
    {
      double result[] = new double[n];
      for (int l = 0; l < n; l++)
      {
        s1[l] = a[l][i];
      } // storing column i in s1

      for (int j = 0; j < i; j++) 
      {
        for (int l2 = 0; l2 < n; l2++)
        {
          s[l2] = q[l2][j];
        }
        double k1 = dot_product(s, s);
        double k2 = dot_product(s, s1);
        for (int i1 = 0; i1 < n; i1++) 
        {
          //result[i1] = result[i1] + ((k2 / (k1 * 1.0)) * s[i1]);
          result[i1] = result[i1] + ((k2*1.0) * s[i1]);
        }

      }
      // storing result of ith column by removing all perpendicular vector which is
      // stored in result matrix
      for (int i2 = 0; i2 < n; i2++)
      {
        q[i2][i] = s1[i2] - result[i2];
      }
     //Dividing matrix by its magnitude 
      double sum=0;
        for (int j1 = 0; j1 < n; j1++)
        {
        sum = sum + (q[j1][i] * q[j1][i]);
        }
      for (int j1 = 0; j1 < n; j1++)
      {
        q[j1][i] = q[j1][i] / Math.sqrt(sum);
      }

    }
  }

  // Function for dot product of two vector.
  public static double dot_product(double a[], double b[]) 
  {
    double sum = 0;
    for (int i = 0; i < a.length; i++) 
    {
      sum = sum + a[i] * b[i];
    }
    return sum;
  }

}