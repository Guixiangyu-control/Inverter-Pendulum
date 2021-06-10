// matrixoperation.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//运行此代码需要安装 Eigen 、Odeint 和 Gnuplot 库。
#include <Eigen/Core>
#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <vector>
#include <cstdio>
#include <boost/numeric/odeint.hpp>
#include <fstream>
#include <iomanip>

#define pi 3.1415926

using namespace std;
using namespace Eigen;
using namespace boost::numeric::odeint;

typedef std::vector<double> state_type; 

void SolveLQR(Matrix<double, 6, 6>& matStateA, Matrix<double, 6, 3>& matStateB,
	Matrix<double, 6, 6>& QofLQR, Matrix<double, 3, 3>& RofLQR, Matrix<double, 6, 6>& PofLQR); //求解LQR
int plotfigure(); //画图函数
//void euqation(const state_type &x, state_type &dx, double t);
void euqation(const state_type &x, state_type &dx, double t); //微分方程
void writetxt();  //为了画图写入txt

struct system  //结构体：系统
{
	Matrix<double, 6, 6 > A;
	Matrix<double, 6, 3 > B;
	Matrix<double, 3, 6 > C;
	Matrix<double, 3, 3 > D;
	Matrix<double, 3, 1> controlledU; /*控制输入*/
} feedbacksystem;

struct state  //结构体：状态变量
{
	state_type theta1;
	state_type theta2;
	state_type theta3;
	state_type theta1_dot;
	state_type theta2_dot;
	state_type theta3_dot;
	state_type time;
} state;


int main()
{
   /* 倒立摆的动力学方程
         J*theta_dot_dot=K*theta+B*Q
   */


	Matrix3d matJ ; //矩阵J

	double J11 = (double)10 / 3;
	double J12 = (double)3 / 2;
	double J13 = (double)1 / 2;
	double J21 = (double)5 / 2;
	double J22 = (double)4 / 3;
	double J23 = (double)1 / 2;
	double J31 = 2;
	double J32 = 1;
	double J33 = (double)7 / 12;
 	
	matJ << J11, J12, J13,
		J21, J22, J23,
		J31, J32, J33;
	 
	Matrix3d matK;  //矩阵K
	double K11 = 25;
	double K12 = 0;
	double K13 = 0;
	double K21 = 0;
	double K22 = 15;
	double K23 = 0;
	double K31 = 0;
	double K32 = 0;
	double K33 = 5;

	matK << K11, K12, K13,
		K21, K22, K23,
		K31, K32, K33;

	Matrix3i matB; //矩阵B
	matB << 1, -1, 0,
		0, 1, -1,
		0, 0, 1;


	//开环系统（倒立摆的）状态空间方程：
	/* 
	         dotx=Ax+Bu
			 y=Cx
	*/
	Matrix<double,6,6> matStateA;   
	Matrix<double,6,3> matStateB;
	Matrix<double,3,6> matStateC;
	

	matStateA << MatrixXd::Zero(3, 3), MatrixXd::Identity(3, 3),
		          matJ.inverse()* matK , MatrixXd::Zero(3, 3);

	
	matStateB << MatrixXd::Zero(3, 3), matJ.inverse()* matK ;


	matStateC << 1, 0, 0, 0, 0, 0,
		         0, 1, 0, 0, 0, 0,
		         0, 0, 1, 0, 0, 0;


	Matrix<double, 6, 6> QofLQR = MatrixXd::Identity(6, 6);  //LQR中的Q

	Matrix<double, 3, 3> RofLQR = MatrixXd::Zero(3, 3);  //LQR中的R

	RofLQR(0, 0) = 10;
	RofLQR(1, 1) = 5;
	RofLQR(2, 2) = 2;
	RofLQR(0, 2) = RofLQR(2, 0) = 1;

	Matrix<double, 6, 6> PofLQR;  //racatti 中的方程中的P
	
	/*racatti方程为 A'P + PA + Q - PBinv(R)B'P=0 */
    SolveLQR(matStateA, matStateB, QofLQR, RofLQR , PofLQR);

	
	Matrix<double, 3, 6> KofLQR; 
    KofLQR = RofLQR.inverse()* matStateB.transpose() * PofLQR; //反馈矩阵

	Matrix<double, 6, 6> feedBackState;
	feedBackState = matStateA - matStateB * KofLQR;  //闭环系统中A为A-B*K

	
	Matrix<double, 3, 3> preMatrix; //前置滤波器
	/*具体可以参考https://zhuanlan.zhihu.com/p/85931427*/
	preMatrix = ( matStateC* (matStateB* KofLQR - matStateA).inverse() * matStateB).inverse();

       
	//反馈系统的状态空间方程
	
	
	feedbacksystem.A = matStateA - matStateB * KofLQR;
	feedbacksystem.B = matStateB * preMatrix;
	feedbacksystem.C = matStateC;
	feedbacksystem.D = MatrixXd::Zero(3,3);

	
	feedbacksystem.controlledU << pi / 2, pi / 2, pi / 2;
	
	const double dt = 0.001;  /*时间间隔*/
	runge_kutta_dopri5<state_type> stepper; /*求解微分方程用runge_kutta算法*/
	state_type x(6);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	x[4] = 0.0;
	x[5] = 0.0;


	

	double t = 0.0;
	state.theta1.push_back(x[0]);
	state.theta2.push_back(x[1]);
	state.theta3.push_back(x[2]);
	state.theta1_dot.push_back(x[3]);
	state.theta2_dot.push_back(x[4]);
	state.theta3_dot.push_back(x[5]);
	state.time.push_back(t);

	for (size_t i(0); i <= 6000; ++i) {
		stepper.do_step( euqation , x, t, dt); //微分方程求解
		t += dt;
		state.theta1.push_back(x[0]);
		state.theta2.push_back(x[1]);
		state.theta3.push_back(x[2]);
		state.theta1_dot.push_back(x[3]);
		state.theta2_dot.push_back(x[4]);
		state.theta3_dot.push_back(x[5]);
		state.time.push_back(t);
	}
	writetxt(); //写入txt
	
	// 关闭文件
	
	plotfigure();
}

/*racatti方程为 A'P + PA + Q - PBinv(R)B'P=0 */
void SolveLQR(Matrix<double, 6, 6>& matStateA, Matrix<double, 6, 3>& matStateB,
	 Matrix<double, 6, 6>& QofLQR, Matrix<double, 3, 3>& RofLQR , Matrix<double, 6, 6>& PofLQR)
{
	const int dim_x = matStateA.rows();
	const int dim_u = matStateB.cols();

	// set Hamilton matrix
	Eigen::MatrixXd Ham = Eigen::MatrixXd::Zero(2 * dim_x, 2 * dim_x);
	Ham << matStateA , -matStateB * RofLQR.inverse() * matStateB.transpose(), -QofLQR, -matStateA.transpose();
	/*  定义汉密尔顿矩阵
	      Ham=[ A  -B*inv(R)*B'
		        -Q  -A']
	*/

	/*参考https://blog.csdn.net/weixin_36815313/article/details/111773535*/
	// calc eigenvalues and eigenvectors
	Eigen::EigenSolver<Eigen::MatrixXd> Eigs(Ham);

	// check eigen values
	// std::cout << "eigen values：\n" << Eigs.eigenvalues() << std::endl;
	// std::cout << "eigen vectors：\n" << Eigs.eigenvectors() << std::endl;

	// extract stable eigenvectors into 'eigvec'
	Eigen::MatrixXcd eigvec = Eigen::MatrixXcd::Zero(2 * dim_x, dim_x);
	int j = 0;
	for (int i = 0; i < 2 * dim_x; ++i) {
		if (Eigs.eigenvalues()[i].real() < 0.) {
			eigvec.col(j) = Eigs.eigenvectors().block(0, i, 2 * dim_x, 1);
			++j;
		}
	}

	// calc P with stable eigen vector matrix
	Eigen::MatrixXcd Vs_1, Vs_2;
	Vs_1 = eigvec.block(0, 0, dim_x, dim_x);
	Vs_2 = eigvec.block(dim_x, 0, dim_x, dim_x);
	PofLQR = (Vs_2 * Vs_1.inverse()).real();
}

int plotfigure()
{
	FILE* pipe = popen("gnuplot", "w"); //具体改为gnuplot的安装位置
	if (pipe == NULL)
	{
		exit(-1);
	}

	fprintf(pipe, "set terminal wxt size 1500, 800\n");
	//fprintf(pipe, "set key fixed left top vertical Right noreverse enhanced autotitle box lt black linewidth 1.000 dashtype solid\n");
	fprintf(pipe, "set multiplot layout 3,2\n");
	//fprintf(pipe, "set title 'Simple Plots'\n");
	//fprintf(pipe, "set title  font ', 20' textcolor lt -1 norotate\n");
	//fprintf(pipe, "set xrange [ * : * ] noreverse writeback\n");
	//fprintf(pipe, "set x2range [ * : * ] noreverse writeback\n");
	//fprintf(pipe, "set yrange [ * : * ] noreverse writeback\n");
	//fprintf(pipe, "set y2range [ * : * ] noreverse writeback\n");
	//fprintf(pipe, "set zrange [ * : * ] noreverse writeback\n");
	//fprintf(pipe, "set cbrange [ * : * ] noreverse writeback\n");
	//fprintf(pipe, "set rrange [ * : * ] noreverse writeback\n");
	//fprintf(pipe, "NO_ANIMATION = 1\n");
	fprintf(pipe, "plot [0:6] [-0.2:1.8] 'state1.d' w l  \n");
	fprintf(pipe, "plot [0:6] [-0.2:1.8] 'state2.d' w l  \n");
	fprintf(pipe, "plot [0:6] [-0.2:1.8] 'state3.d' w l  \n");
	fprintf(pipe, "plot [0:6] [-0.6:1.8] 'state4.d' w l  \n");
	fprintf(pipe, "plot [0:6] [-0.2:1.8] 'state5.d' w l  \n");
	fprintf(pipe, "plot [0:6] [-0.2:1.8] 'state6.d' w l  \n");

	fprintf(pipe, "pause mouse\n");
	
        pclose(pipe);


	

	return 0;
}




void euqation(const state_type &x , state_type &dx , double t ) //系统微分方程
{
	
	dx[0] = feedbacksystem.A(0, 0)*x[0] + feedbacksystem.A(0, 1)*x[1] + feedbacksystem.A(0, 2)*x[2] + feedbacksystem.A(0, 3)*x[3] 
		+ feedbacksystem.A(0, 4)*x[4] + feedbacksystem.A(0, 5)* x[5] + feedbacksystem.B(0, 0)*feedbacksystem.controlledU(0, 0) 
		+ feedbacksystem.B(0, 1)*feedbacksystem.controlledU(1, 0) + feedbacksystem.B(0, 2)*feedbacksystem.controlledU(2, 0);

	dx[1] = feedbacksystem.A(1, 0)*x[0] + feedbacksystem.A(1, 1)*x[1] + feedbacksystem.A(1, 2)*x[2] + feedbacksystem.A(1, 3)*x[3] 
		+ feedbacksystem.A(1, 4)*x[4] + feedbacksystem.A(1, 5)* x[5] + feedbacksystem.B(1, 0)*feedbacksystem.controlledU(0, 0) 
		+ feedbacksystem.B(1, 1)*feedbacksystem.controlledU(1, 0) + feedbacksystem.B(1, 2)*feedbacksystem.controlledU(2, 0);

	dx[2] = feedbacksystem.A(2, 0)*x[0] + feedbacksystem.A(2, 1)*x[1] + feedbacksystem.A(2, 2)*x[2] + feedbacksystem.A(2, 3)*x[3] 
		+ feedbacksystem.A(2, 4)*x[4] + feedbacksystem.A(2, 5)* x[5] + feedbacksystem.B(2, 0)*feedbacksystem.controlledU(0, 0) 
		+ feedbacksystem.B(2, 1)*feedbacksystem.controlledU(1, 0) + feedbacksystem.B(2, 2)*feedbacksystem.controlledU(2, 0);

	dx[3] = feedbacksystem.A(3, 0)*x[0] + feedbacksystem.A(3, 1)*x[1] + feedbacksystem.A(3, 2)*x[2] + feedbacksystem.A(3, 3)*x[3] 
		+ feedbacksystem.A(3, 4)*x[4] + feedbacksystem.A(3, 5)* x[5] + feedbacksystem.B(3, 0)*feedbacksystem.controlledU(0, 0) 
		+ feedbacksystem.B(3, 1)*feedbacksystem.controlledU(1, 0) + feedbacksystem.B(3, 2)*feedbacksystem.controlledU(2, 0);

	dx[4] = feedbacksystem.A(4, 0)*x[0] + feedbacksystem.A(4, 1)*x[1] + feedbacksystem.A(4, 2)*x[2] + feedbacksystem.A(4, 3)*x[3] 
		+ feedbacksystem.A(4, 4)*x[4] + feedbacksystem.A(4, 5)* x[5] + feedbacksystem.B(4, 0)*feedbacksystem.controlledU(0, 0) 
		+ feedbacksystem.B(4, 1)*feedbacksystem.controlledU(1, 0) + feedbacksystem.B(4, 2)*feedbacksystem.controlledU(2, 0);

	dx[5] = feedbacksystem.A(5, 0)*x[0] + feedbacksystem.A(5, 1)*x[1] + feedbacksystem.A(5, 2)*x[2] + feedbacksystem.A(5, 3)*x[3] 
		+ feedbacksystem.A(5, 4)*x[4] + feedbacksystem.A(5, 5)* x[5] + feedbacksystem.B(5, 0)*feedbacksystem.controlledU(0, 0) 
		+ feedbacksystem.B(5, 1)*feedbacksystem.controlledU(1, 0) + feedbacksystem.B(5, 2)*feedbacksystem.controlledU(2, 0);

}

void writetxt()
{
	ofstream out_txt_file;
	
	//写入状态1
	out_txt_file.open("./state1.d", ios::out | ios::trunc);
	out_txt_file << "#" << "X" << " " << " " << "Y" << endl;
	for (int i = 0; i < state.theta1.size(); i++)
	{
		out_txt_file << "  " << setprecision(4) << state.time[i] << " " << setprecision(4) << state.theta1[i] << endl;
	}
	out_txt_file.close();

	ofstream out_txt_file1;
	//写入状态2
	out_txt_file1.open("./state2.d", ios::out | ios::trunc);
	out_txt_file1 << "#" << "X" << " " << " " << "Y" << endl;
	for (int i = 0; i < state.theta2.size(); ++i)
	{
		out_txt_file1 << "  " << setprecision(4) << state.time[i] << " " << setprecision(4) << state.theta2[i] << endl;
	}
	out_txt_file1.close();

	ofstream out_txt_file2;
	//写入状态3
	out_txt_file2.open("./state3.d", ios::out | ios::trunc);
	out_txt_file2 << "#" << "X" << " " << " " << "Y" << endl;
	for (int i = 0; i < state.theta3.size(); ++i)
	{
		out_txt_file2 << "  " << setprecision(4) << state.time[i] << " " << setprecision(4) << state.theta3[i] << endl;
	}
	out_txt_file2.close();

	ofstream out_txt_file3;
	//写入状态4
	out_txt_file3.open("./state4.d", ios::out | ios::trunc);
	out_txt_file3 << "#" << "X" << " " << " " << "Y" << endl;
	for (int i = 0; i < state.theta1_dot.size(); i++)
	{
		out_txt_file3 << "  " << setprecision(4) << state.time[i] << " " << setprecision(4) << state.theta1_dot[i] << endl;
	}
	out_txt_file3.close();

	ofstream out_txt_file4;
	//写入状态5
	out_txt_file4.open("./state5.d", ios::out | ios::trunc);
	out_txt_file4 << "#" << "X" << " " << " " << "Y" << endl;
	for (int i = 0; i < state.theta2_dot.size(); i++)
	{
		out_txt_file4 << "  " << setprecision(4) << state.time[i] << " " << setprecision(4) << state.theta2_dot[i] << endl;
	}
	out_txt_file4.close();


	ofstream out_txt_file5;
	//写入状态6
	out_txt_file5.open("./state6.d", ios::out | ios::trunc);
	out_txt_file5 << "#" << "X" << " " << " " << "Y" << endl;
	for (int i = 0; i < state.theta3_dot.size(); i++)
	{
		out_txt_file5 << "  " << setprecision(4) << state.time[i] << " " << setprecision(4) << state.theta3_dot[i] << endl;
	}

	out_txt_file5.close();


}
