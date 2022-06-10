#include "ScanPathGeneration.h"
#include<math.h>
#include<time.h>
using namespace MeshLib;

double pointDistance(pathPoint p1, pathPoint p2)
{
	return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
}
double normalAngle(pathPoint p1, pathPoint p2)
{
	double cosA;
	cosA = (p1.nx*p2.nx + p1.ny*p2.ny + p1.nz*p2.nz) / (sqrt(p1.nx*p1.nx + p1.ny*p1.ny + p1.nz*p1.nz)*sqrt(p2.nx*p2.nx + p2.ny*p2.ny + p2.nz*p2.nz));
	return cosA;
}
Point normal(pathPoint p1, pathPoint p2)
{
	double x,y,z,s;
	s = sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
	x = (p2.x - p1.x) / s;
	y = (p2.y - p1.y) / s;
	z = (p2.z - p1.z) / s;
	Point n(x, y, z);
	return n;
}
double normalAngle(Point p1, Point p2)
{
	double cosA;
	cosA = (p1.x()*p2.x() + p1.y()*p2.y() + p1.z()*p2.z()) / (sqrt(p1.x()*p1.x() + p1.y()*p1.y() + p1.z()*p1.z())*sqrt(p2.x()*p2.x() + p2.y()*p2.y() + p2.z()*p2.z()));
	return cosA;
}
SPGen::SPGen(Mesh *mesh, double s)
{
	scanInterval = s;
	for (MeshVertexIterator viter(mesh); !viter.end(); ++viter) {
		Vertex * v = *viter;
		path_verts.push_back(v);
	}
	verts2points();

	bool isFinished=false;
	pathPoint p;
	double dist, cosA;
	int k=3;
	//while (!isFinished)//sampling 
	//{
	//	isFinished = true;
	//	for (int i = 0; i < pathTemp.size(); i++)
	//	{
	//		p = pathTemp[i];
	//		if (!p.is_deleted)
	//		{
	//			//if (p.is_boundary)
	//			if (0)
	//			{
	//				for (int j = 0; j < pathTemp.size(); j++)
	//				{
	//					if (!pathTemp[j].is_deleted&&j!=i)
	//					{
	//						dist = pointDistance(p, pathTemp[j]);
	//						cosA = normalAngle(p, pathTemp[j]);
	//						if (cosA >= 0) cosA = pow(cosA, k);
	//						else cosA = 0;
	//						if (dist < scanInterval*cosA)
	//						{
	//							//if (cosA < 0.7)
	//							//{
	//							//	std::cout << dist << "   " << cosA << std::endl;
	//							//}
	//							
	//							pathTemp[j].is_deleted = true;
	//							isFinished = false;
	//						}

	//					}

	//				}
	//			}
	//			else
	//			{
	//				for (int j = 0; j < pathTemp.size(); j++)
	//				{
	//					//if (pathTemp[j].is_boundary)
	//						;
	//					if (0);
	//					else
	//					{
	//						if (!pathTemp[j].is_deleted&&j != i)
	//						{
	//							dist = pointDistance(p, pathTemp[j]);
	//							//cosA = pow(normalAngle(p, scanPath[j]), 2 * k);
	//							cosA = normalAngle(p, pathTemp[j]);
	//							if (cosA >= 0) cosA = pow(cosA, k);
	//							else cosA = 0;
	//							if (dist < scanInterval*cosA)
	//							{
	//								//if (cosA < 0.7)
	//								//{
	//								//	std::cout << dist << "   " << cosA << std::endl;
	//								//}

	//								pathTemp[j].is_deleted = true;
	//								isFinished = false;
	//							}

	//						}
	//					}

	//				}
	//			}
	//		}
	//	}
	//}
	// 行扫法采样
	double dd,scanAxis;
	double interval = 50;//曲线中采样点间距
	for (int i = 0; i < pathTemp.size(); i++)
	{
		scanAxis = pathTemp[i].z;
		dd = fmod(scanAxis,scanInterval);
		if (abs(dd) < 10) ;
		else
			pathTemp[i].is_deleted = true;
	}
	while (!isFinished)//sampling 
	{
		isFinished = true;
		for (int i = 0; i < pathTemp.size(); i++)
		{
			p = pathTemp[i];
			if (!p.is_deleted)
			{
				//if (p.is_boundary)
				if (0)
				{
					for (int j = 0; j < pathTemp.size(); j++)
					{
						if (!pathTemp[j].is_deleted&&j != i)
						{
							dist = pointDistance(p, pathTemp[j]);
							cosA = normalAngle(p, pathTemp[j]);
							if (cosA >= 0) cosA = pow(cosA, k);
							else cosA = 0;
							if (dist < interval*cosA)
							{
								//if (cosA < 0.7)
								//{
								//	std::cout << dist << "   " << cosA << std::endl;
								//}

								pathTemp[j].is_deleted = true;
								isFinished = false;
							}

						}

					}
				}
				else
				{
					for (int j = 0; j < pathTemp.size(); j++)
					{
						//if (pathTemp[j].is_boundary)
						;
						if (0);
						else
						{
							if (!pathTemp[j].is_deleted&&j != i)
							{
								dist = pointDistance(p, pathTemp[j]);
								//cosA = pow(normalAngle(p, scanPath[j]), 2 * k);
								cosA = normalAngle(p, pathTemp[j]);
								if (cosA >= 0) cosA = pow(cosA, k);
								else cosA = 0;
								if (dist < interval*cosA)
								{
									//if (cosA < 0.7)
									//{
									//	std::cout << dist << "   " << cosA << std::endl;
									//}

									pathTemp[j].is_deleted = true;
									isFinished = false;
								}

							}
						}

					}
				}
			}
			
		}
	}
	std::cout << "--> Writing sampled vertexes..." << std::endl;
	writerVertexes("../output/0_vertexes.txt");//输出采样点

 	//openFile("../output/0_vertexes_totle.txt");//输入采样点进行后续操作

	offset();
	
	stationPlan();
	genBox(1100,500);

	writer3dPath("../output/1_original_path.txt");

	//writerPath("../output/2_projection_path.txt");
	//writePolyline("../output/3_original_box.poly");
	//filtBox(120);
	//writePolyline("../output/4_sampled_box.poly");
	////filter();

	//genStation();
	//writeStationBox("../output/5_stationsBox.poly");
	//writeStation("../output/6_stations.txt");
	//writeStationPath("../output/7_stationsPath.txt");
	//inorder(2);
 //   writepoly("../output/8_path_s1.poly");
}

SPGen::~SPGen(){}
void SPGen::openFile(const char *filename)
{
	pathTemp.clear();//清空
	pathPoint p;
	FILE *stream;//open the file and take the data
	stream = fopen(filename, "r");
	int P_num = 0;
	while (!feof(stream))
	{

		if (9 == fscanf(stream, "%lf%lf%lf%lf%lf%lf%lf%lf%lf", &p.x, &p.y, &p.z, &p.r, &p.g, &p.b, &p.nx, &p.ny, &p.nz))
		{
			pathTemp.push_back(p);
			P_num++;
		}
	}
	std::cout << "Data read finished." << std::endl;
	std::cout << std::endl;
	std::cout << "Totle:" << P_num << std::endl;
	fclose(stream);//close the file
}
void SPGen::writerVertexes(const char *output)
{
	FILE *stream;//open the file and take the data
	errno_t err;
	err = fopen_s(&stream, output, "w");
	int i = 0, num = 0;
	if (err == 0)
	{
		for (; i < pathTemp.size(); i++)
		{
			if (!pathTemp[i].is_deleted)
			{
				//fprintf_s(stream, "%lf %lf %lf %lf %lf %lf %i %i %i\n", pathTemp[i].x, pathTemp[i].y, pathTemp[i].z, pathTemp[i].nx, pathTemp[i].ny, pathTemp[i].nz, pathTemp[i].r, pathTemp[i].g, pathTemp[i].b);
				fprintf_s(stream, "%lf %lf %lf %lf %lf %lf\n", pathTemp[i].x, pathTemp[i].y, pathTemp[i].z, pathTemp[i].nx, pathTemp[i].ny, pathTemp[i].nz);
				num += 1;
			}

		}
		err = fclose(stream);//close the file
		if (err == 0)
			std::cout << "Total path points:" << num << std::endl << "Saved sucessfully!\n";
		else
			std::cout << "erro: file could not closed\n";
	}
	else
		std::cout << "erro: file could not opened\n";
}
void SPGen::writerPath(const char *output)
{
	FILE *stream;//open the file and take the data
	errno_t err;
	err = fopen_s(&stream, output, "w");
	int i = 0,num=0;
	if (err == 0)
	{
		for (; i < scanPath.size(); i++)
		{
			if (!scanPath[i].is_deleted)
			{
				//fprintf_s(stream, "%lf %lf %lf %lf %lf %lf %i %i %i\n", scanPath[i].x, scanPath[i].y, scanPath[i].z, scanPath[i].nx, scanPath[i].ny, scanPath[i].nz, scanPath[i].r, scanPath[i].g, scanPath[i].b);
				fprintf_s(stream, "%lf %lf %lf %lf %lf %lf\n", scanPath[i].x, scanPath[i].y, scanPath[i].z, scanPath[i].nx, scanPath[i].ny, scanPath[i].nz);
				num += 1;
			}
			    
		}
		err = fclose(stream);//close the file
		if (err == 0)
			std::cout << "Total path points:"<<num<<std::endl<<"Saved sucessfully!\n";
		else
			std::cout << "erro: file could not closed\n";
	}
	else
		std::cout << "erro: file could not opened\n";
}
void SPGen::writer3dPath(const char *output)
{
	FILE *stream;//open the file and take the data
	errno_t err;
	err = fopen_s(&stream, output, "w");
	int i = 0, num = 0;
	if (err == 0)
	{
		for (; i < outPath.size(); i++)
		{
			if (!scanPath[i].is_deleted)
			{
				//fprintf_s(stream, "%lf %lf %lf %lf %lf %lf %i %i %i\n", scanPath[i].x, scanPath[i].y, scanPath[i].z, scanPath[i].nx, scanPath[i].ny, scanPath[i].nz, scanPath[i].r, scanPath[i].g, scanPath[i].b);
				fprintf_s(stream, "%lf %lf %lf %lf %lf %lf\n", outPath[i].x, outPath[i].y, outPath[i].z, outPath[i].nx, outPath[i].ny, outPath[i].nz);
				num += 1;
			}

		}
		err = fclose(stream);//close the file
		if (err == 0)
			std::cout << "Total path points:" << num << std::endl << "Saved sucessfully!\n";
		else
			std::cout << "erro: file could not closed\n";
	}
	else
		std::cout << "erro: file could not opened\n";
}
void SPGen::writepoly(const char *filename)
{
	FILE *stream;//open the file and take the data
	errno_t err;

	err = fopen_s(&stream,filename, "w");
	int i = 0, num = 0;
	if (err == 0)
	{

		for (; i < outPath.size(); i++)
		{
			if (!outPath[i].is_deleted)
			{
				//fprintf_s(stream, "%lf %lf %lf %lf %lf %lf %i %i %i\n", scanPath[i].x, scanPath[i].y, scanPath[i].z, scanPath[i].nx, scanPath[i].ny, scanPath[i].nz, scanPath[i].r, scanPath[i].g, scanPath[i].b);
				fprintf_s(stream, "%lf %lf %lf\n", outPath[i].x, outPath[i].y, outPath[i].z);
				num += 1;
			}


		}
		err = fclose(stream);//close the file
		if (err == 0)
			std::cout << "Total path points:" << num << std::endl << "Saved sucessfully!\n";
		else
			std::cout << "erro: file could not closed\n";
	}
	else
		std::cout << "erro: file could not opened\n";
}
void SPGen::writePolyline(const char *filename)
{
	FILE *stream;//open the file and take the data
	errno_t err;
	std::vector<pathPoint> temp;
	err = fopen_s(&stream, filename, "w");
	int num = 0;
	if (err == 0)
	{

		for (int i = 0; i < polyLine.size(); i++)
		{
			temp = polyLine[i];
			if (!temp[0].is_deleted)
			{
				for (int j = 0; j < temp.size(); j++)
				{
					if (j != 2)
						fprintf_s(stream, "%lf %lf %lf\n", temp[j].x, temp[j].y, temp[j].z);
				}
				fprintf_s(stream, "\n");
				num += 1;
			}
		}
		err = fclose(stream);//close the file
		if (err == 0)
			std::cout << "Total polylines:" << num << std::endl <<" -> "<<filename<< "\n";
		else
			std::cout << "erro: file could not closed\n";
	}
	else
		std::cout << "erro: file could not opened\n";
}
void SPGen::writeStationBox(const char *filename)
{
	FILE *stream;//open the file and take the data
	errno_t err;
	std::vector<pathPoint> temp;
	err = fopen_s(&stream, filename, "w");
	int num = 0;
	if (err == 0)
	{

		for (int i = 0; i < station.size(); i++)
		{
			temp = polyLine[station[i]];
			if (!temp[0].is_deleted)
			{
				for (int j = 0; j < temp.size(); j++)
				{
					if(j!=2)
						fprintf_s(stream, "%lf %lf %lf\n", temp[j].x, temp[j].y, temp[j].z);
				}
				fprintf_s(stream, "\n");
				num += 1;
			}
		}
		err = fclose(stream);//close the file
		if (err == 0)
			std::cout << "Total polylines:" << num << std::endl << " -> " << filename << "\n";
		else
			std::cout << "erro: file could not closed\n";
	}
	else
		std::cout << "erro: file could not opened\n";
}

void SPGen::writeStation(const char *filename)
{
	FILE *stream;//open the file and take the data
	errno_t err;
	std::vector<pathPoint> temp;
	err = fopen_s(&stream, filename, "w");
	int num = 0;
	double set = 500;
	if (err == 0)
	{

		for (int i = 0; i < station.size(); i++)
		{
			std::cout << "s: " << station[i] << std::endl;
			temp = polyLine[station[i]];
			if (!temp[0].is_deleted)
			{
				num += 1;
				fprintf_s(stream, "%lf %lf %lf\n", temp[2].x - set * temp[2].nx, temp[2].y, temp[2].z - set * temp[2].nz);
			}
		}
		err = fclose(stream);//close the file
		if (err == 0)
			std::cout << "Total polylines:" << num << std::endl << " -> " << filename << "\n";
		else
			std::cout << "erro: file could not closed\n";
	}
	else
		std::cout << "erro: file could not opened\n";
}
void SPGen::writeStationPath(const char *filename)
{
	FILE *stream;//open the file and take the data
	errno_t err;
	std::vector<int> pointList;
	pathPoint temp;
	err = fopen_s(&stream, filename, "w");
	int num = 0;
	if (err == 0)
	{

		for (int i = 0; i < station.size(); i++)
		{
			fprintf_s(stream, "%i\n",station[i]);
			pointList = coverPoint[station[i]];
			for (int m = 0; m < pointList.size(); m++)
			{
				temp = outPath[pointList[m]];
				if (!temp.is_deleted)
				{
					fprintf_s(stream, "%lf %lf %lf %lf %lf %lf\n", temp.x, temp.y, temp.z, temp.nx, temp.ny, temp.nz);
					num += 1;
				}
			}
			
		}
		err = fclose(stream);//close the file
		if (err == 0)
			std::cout << "Total polylines:" << num << std::endl << " -> " << filename << "\n";
		else
			std::cout << "erro: file could not closed\n";
	}
	else
		std::cout << "erro: file could not opened\n";
}

void SPGen::verts2points()
{
	int num = 0;
	for (std::list<Vertex*>::iterator viter = path_verts.begin(); viter != path_verts.end(); ++viter)
	{
		Vertex     *v = *viter;
		pathPoint path;
		Point p, n;
		bool b;
		p = v->point();
		n = v->normal();
		b = v->boundary();
		path.x = p.x();
		path.y = p.y();
		path.z = p.z();
		path.nx = n.x();
		path.ny = n.y();
		path.nz = n.z();
		path.is_boundary = b;
		pathTemp.push_back(path);
		num += 1;
	}
	std::cout << "original vertics points total:" << num << std::endl;
}
void SPGen::offset()//采样点沿法线方向偏置得到路径点
{
	double s;
	double xc=0, yc=0, zc=0, nxc=0, nyc=0, nzc=0;
	int num = 0;
	pathPoint p;
	for (int i=0; i < pathTemp.size(); i++)
	{
		if (!pathTemp[i].is_deleted)
		{

			s = sqrt(pathTemp[i].nx*pathTemp[i].nx + pathTemp[i].ny*pathTemp[i].ny + pathTemp[i].nz*pathTemp[i].nz);
			//p.x = pathTemp[i].x + offsetDist * pathTemp[i].nx / s;
			//p.y = pathTemp[i].y + offsetDist * pathTemp[i].ny / s;
			//p.z = pathTemp[i].z + offsetDist * pathTemp[i].nz / s;
			//p.nx = pathTemp[i].nx * -1;
			//p.ny = pathTemp[i].ny * -1;
			//p.nz = pathTemp[i].nz * -1;
			p.x = pathTemp[i].x - offsetDist * pathTemp[i].nx / s;
			p.y = pathTemp[i].y - offsetDist * pathTemp[i].ny / s;
			p.z = pathTemp[i].z - offsetDist * pathTemp[i].nz / s;
			p.nx = pathTemp[i].nx * 1;
			p.ny = pathTemp[i].ny * 1;
			p.nz = pathTemp[i].nz * 1;
			xc += p.x;//计算路径点几何中心
			yc += p.y;
			zc += p.z;
			nxc += p.nx;
			nyc += p.ny;
			nzc += p.nz;
			num++;
			scanPath.push_back(p);
			pathPoint_left.push_back(true);
		}

	}
	std::cout << "scan path points total:" << num << std::endl;
	path_center.x = xc / num;//计算路径点几何中心
	path_center.y = yc / num;
	path_center.z = zc / num;
	path_center.nx = nxc / num;
	path_center.ny = nyc / num;
	path_center.nz = nzc / num;
}
void SPGen::inorder(int w)
{
	double cos_alpha, cos_beta;
	int pm = 0.95;
	double k = 0.02;
	double d,d0 = 3*scanInterval;
	double s,max_s;
	int nextPoint;
	pathPoint p0, p1;
	Point n0, n1, initial_n(0, 0, 1);
	std::vector<pathPoint> temp;
	bool finished = false;
	//int color = 30;
	double d_s, a_s, b_s;
	std::vector<int> pointList;
	pathPoint t;
	std::vector<pathPoint> tempPath;

	tempPath.clear();
	pointList = coverPoint[station[w]];
	for (int m = 0; m < pointList.size(); m++)
	{
		t = outPath[pointList[m]];
		if (!t.is_deleted)
		{
			tempPath.push_back(t);
		}
	}

	//scanPath[0].r = 255;
	//scanPath[0].g = 255;
	//scanPath[0].b = 255;
	tempPath[0].is_inorder = true;
	temp.push_back(tempPath[0]);
	n0 = initial_n;
	p0 = tempPath[0];
	while (!finished)
	{
		finished = true;
		for (int i = 1; i < tempPath.size(); i++)
		{
			if (!tempPath[i].is_deleted)
			{
				if (!tempPath[i].is_inorder)
				{
					p1 = tempPath[i];
					n1 = normal(p0, p1);
					cos_alpha = normalAngle(n1, n0);
					cos_beta = normalAngle(p1, p0);
					d = pointDistance(p1, p0);
					//s = (pm * (1 + cos_alpha) + (1 - pm)*(1 + cos_beta)) / (2 * (1 + exp(k*(d - d0))));
					d_s = 1 / (1 + exp(k*(d - d0)));
					//a_s = (1 + pm * cos_alpha);
					//b_s = (1 + (1 - pm)*cos_beta);
					//s = a_s *b_s*d_s;
					a_s = pm * (1 + cos_alpha);
					b_s = (1 - pm)*(1 + cos_beta);
					s = (a_s +b_s)*d_s;
					if (finished|| s >max_s)
					{
 						max_s = s;
						nextPoint = i;
					}
					finished = false;
				}
			}
		}

		//tempPath[nextPoint].r = color % 255;
		//color += 30;
		tempPath[nextPoint].is_inorder = true;
		temp.push_back(tempPath[nextPoint]);
		n0 = normal(p0, tempPath[nextPoint]);
		p0 = temp.back();
	}

	outPath.clear();

	for (int j=0; j < temp.size(); j++)
	{
		outPath.push_back(temp[j]);
	}
}
void SPGen::stationPlan()//didn't finished
{
	//std::vector<std::vector<pathPoint>> S;
	//std::vector<pathPoint> s0,s1, s2;
	//int n = 0;
	double s;
	for (int i = 0; i < scanPath.size(); i++)
	{
		if (!scanPath[i].is_deleted)
		{
			outPath.push_back(scanPath[i]);
			scanPath[i].y = 0;//将y轴压缩，即投影至xz平面
			scanPath[i].ny = 0;
			s = sqrt(scanPath[i].nx*scanPath[i].nx + scanPath[i].nz*scanPath[i].nz);
			scanPath[i].nx = scanPath[i].nx / s;
			scanPath[i].nz = scanPath[i].nz / s;
		}
	}
}
void SPGen::genBox(double l,double w)
{
	std::vector<pathPoint> box;
	std::vector<int> pointList;
	pathPoint p1, p2, p3, p4,p5;
	double tanTheta;
	double dx1, dz1,dx2,dz2;
	double d;
	int num=0;
	bool inbox;
	
	for (int i = 0; i < scanPath.size(); i++)
	{
		if (!scanPath[i].is_deleted)
		{
			//判断其余点是否在框内
			pointList.clear();
			for (int k = 0; k < scanPath.size(); k++)
			{
				if (!scanPath[k].is_deleted)
				{
					inbox = is_inbox(scanPath[i], scanPath[k], l, w);
					if (inbox)
					{
						pointList.push_back(k);
						scanPath[i].cover++;
					}
				}
			}
			if (scanPath[i].cover < 10||normalAngle(scanPath[i],path_center)<0)//删除异常点
			{
				scanPath[i].is_deleted = true;
				pathPoint_left[i] = false;
				
				continue;
			}
			box.clear();
			//绘制框的四个点，连成一个矩形框
			tanTheta = scanPath[i].nx / scanPath[i].nz;
			dx1 = (l / 2 + w / 2 * tanTheta)*scanPath[i].nz;
			dz1 = w / 2 / scanPath[i].nz - dx1 * tanTheta;
			dx2 = (l / 2 - w / 2 * tanTheta)*scanPath[i].nz;
			dz2= w / 2 / scanPath[i].nz + dx2 * tanTheta;
			p1.x = scanPath[i].x+dx1;
			p1.y = scanPath[i].y;
			p1.z = scanPath[i].z+dz1;
			p1.cover = scanPath[i].cover;
			box.push_back(p1);
			p2.x = scanPath[i].x + dx2;
			p2.y = scanPath[i].y;
			p2.z = scanPath[i].z - dz2;
			box.push_back(p2);
			p3.x = scanPath[i].x;
			p3.y = scanPath[i].y;
			p3.z = scanPath[i].z;
			p3.nx = scanPath[i].nx;
			p3.ny = scanPath[i].ny;
			p3.nz = scanPath[i].nz;
			box.push_back(p3);
			p4.x = scanPath[i].x - dx1;
			p4.y = scanPath[i].y;
			p4.z = scanPath[i].z - dz1;
			box.push_back(p4);
			p5.x = scanPath[i].x - dx2;
			p5.y = scanPath[i].y;
			p5.z = scanPath[i].z + dz2;
			box.push_back(p5);
			box.push_back(p1);

			std::cout << "polyline"<< num<<" cover:" << scanPath[i].cover << std::endl;
			polyLine.push_back(box);
			box_left.push_back(true);
			coverPoint.push_back(pointList);
			
			
			num++;
		}
	}
	std::cout << "num:"<<num << std::endl;
}
void SPGen::filtBox(double d)
{
	std::vector<pathPoint> box1,box2;
	for (int i = 0; i < polyLine.size(); i++)
	{
		box1 = polyLine[i];
		if (!box1[0].is_deleted)
		{
			for (int j = 0; j < polyLine.size(); j++)
			{
				if (i != j)
				{
					box2 = polyLine[j];
					if (!box2[0].is_deleted)
					{
						if (pointDistance(box1[3], box2[3]) < d)
						{
							if (box1[0].cover >= box2[0].cover)
							{
								box2[0].is_deleted = true;
								polyLine[j] = box2;
								box_left[j] = false;
							}
							else
							{
 								box1[0].is_deleted = true;
								polyLine[i] = box1;
								box_left[j] = false;
								continue;
							}
						}
					}
				}


			}
		}
		
	}
}
void SPGen::genStation()
{
	double nc = 0;
	double r1, r2;
	double d;
	double theta;
	double s;
	double s_min=10000;
	clock_t start, finish;
	double duration;

	std::vector<pathPoint> box1, box2,box3;
	station.push_back(0);
	station.push_back(0);
	station.push_back(0);
	
	for (int i = 0; i < polyLine.size(); i++)
	{
		box1 = polyLine[i];
		if (!box1[0].is_deleted)
		{
			start = clock();
			for (int j = 0; j < polyLine.size(); j++)
			{
				if (i != j)
				{
					box2 = polyLine[j];
					if (!box2[0].is_deleted)
					{
						for (int k = 0; k < polyLine.size(); k++)
						{
							if (k != j && k != i)
							{
								box3 = polyLine[k];
								if (!box3[0].is_deleted)
								{
									nc = numOfNoCovered(i, j, k);
									r1 = pointDistance(box1[2], box2[2])/1000;
									r2 = pointDistance(box2[2], box3[2])/1000;
									//d = double(1 / r1) + double(1 / r2) ;
									d = abs(r1-1.1) + abs(r2-1.1);
									theta =2- normalAngle(box1[2], box2[2])- normalAngle(box2[2], box3[2]);
									s = 1*nc + 1*d+theta;
									//s = 1 * nc + 1 * d;
									if (s < s_min)
									{
										s_min = s;
										std::cout << std::endl;
										std::cout << " i : " << i << "  j : " << j << "  k: " << k << std::endl;
										std::cout << " s_min =" << s_min  << std::endl;

										station[0]=i;
										station[1] = j;
										station[2] = k;
									}
								}
							}
						}
					}
				}
			}
			finish = clock();
			duration = double(finish - start) / CLOCKS_PER_SEC;
			std::cout << " *** " << i << " : " << duration << std::endl;
		}
	}

}
int SPGen::numOfNoCovered(int i,int j,int k)
{
	int num=0;
	std::vector<int> pointList;
	bool in;
	for (int l = 0; l < pathPoint_left.size(); l++)
	{
		in = false;
		if (pathPoint_left[l])
		{
			if (!in)
			{
				pointList = coverPoint[i];
				for (int m = 0; m < pointList.size(); m++)
				{
					if (pointList[m] == l)
					{
						in = true;
						continue;
					}
				}
			}
			if (!in)
			{
				pointList = coverPoint[j];
				for (int m = 0; m < pointList.size(); m++)
				{
					if (pointList[m] == l)
					{
						in = true;
						continue;
					}
				}
			}
			if (!in)
			{
				pointList = coverPoint[k];
				for (int m = 0; m < pointList.size(); m++)
				{
					if (pointList[m] == l)
					{
						in = true;
						continue;
					}
				}
			}
			if (!in) num++;
		}
	}
	return num;
}
bool SPGen::is_inbox(pathPoint box_center, pathPoint p, double l, double w)
{
	double x0 = box_center.x;
	double y0 = box_center.z;
	double cosT = box_center.nz;
	double sinT = box_center.nx;
	Eigen::Vector3d pc, p0, p1;
	pc << box_center.x, box_center.z, 1;
	p0 << p.x, p.z, 1;
	double x_max = pc[0] + l / 2;
	double x_min = pc[0] - l / 2;
	double y_max = pc[1] + w / 2;
	double y_min = pc[1] - w / 2;

	
	Eigen::Matrix3d t1,t2,R,M;
	t1 << 1, 0, x0, 0, 1, y0, 0, 0, 1;
	t2 << 1, 0, -x0, 0, 1, -y0, 0, 0, 1;
	R << cosT, -sinT, 0, sinT, cosT, 0, 0, 0, 1;
	M = t1 * R*t2;
	p1 = M * p0;
	double cosn = normalAngle(box_center, p);

	if (p1[0] > x_min&&p1[0] < x_max&&p1[1] < y_max&&p1[1] > y_min&&cosn>0)
		return true;
	else
		return false;
}
double SPGen::IoU(std::vector<pathPoint> box1, std::vector<pathPoint> box2)
{
	return 0;
}
void SPGen::filter()
{
	double a = 0;
	double b = 2600;
	for (int i = 0; i < scanPath.size(); i++)
	{
		if (scanPath[i].z<a|| scanPath[i].z > b)
		{
			scanPath[i].is_deleted = true;
		}

	}
}
//std::vector<pathPoint> SPGen::verts2points(std::list<Vertex*> pv)
//{
//	std::vector<pathPoint> sPath;
//	for (std::list<Vertex*>::iterator viter = pv.begin(); viter != pv.end(); ++viter)
//	{
//		Vertex     *v = *viter;
//		pathPoint path;
//		Point p, n;
//		bool b;
//		p = v->point();
//		n = v->normal();
//		b = v->boundary();
//		path.x =p.x();
//		path.y = p.y();
//		path.z = p.z();
//		path.nx = n.x();
//		path.ny = n.y();
//		path.nz = n.z();
//		path.is_boundary = b;
//		sPath.push_back(path);
//	}
//	return sPath;
//}
pathPoint SPGen::vert2point(Vertex* v)
{
	pathPoint path;
	Point p, n;
	bool b;
	p = v->point();
	n = v->normal();
	b = v->boundary();
	path.x = p.x();
	path.y = p.y();
	path.z = p.z();
	path.nx = n.x();
	path.ny = n.y();
	path.nz = n.z();
	path.is_boundary = b;
	return path;
}

//double pointDistance(pathPoint p1, pathPoint p2)
//{
//	return sqrt((p1.x - p2.x)*(p1.x - p2.x) + (p1.y - p2.y)*(p1.y - p2.y) + (p1.z - p2.z)*(p1.z - p2.z));
//}
//
//double normalAngle(pathPoint p1, pathPoint p2)
//{
//	return abs((p1.nx*p2.nx + p1.ny*p2.ny + p1.nz*p2.nz) / (sqrt(p1.nx*p1.nx + p1.ny*p1.ny + p1.nz*p1.nz)*sqrt(p2.nx*p2.nx + p2.ny*p2.ny + p2.nz*p2.nz)));
//}