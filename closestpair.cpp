#include <iostream>
#include <string>
#include <cstdlib>
#include <vector>
#include <tuple>
#include <cmath>
#include <time.h>
#include <math.h>
#include <chrono>

using namespace std;

//naive way to get to get closestPair of points
void closest_naive(const char* input)
{

	//array of Points each with the second part of the array being x and y values
	int points [stoi(input)][2];

	//populate the points with random values between o and a milion
	for(int k = 0; k < stoi(input); k++)
	{
		points[k][0] = rand() %1000000;
		points[k][1] = rand() %1000000;
	}

	//print out the pairs of values
	//for(int k = 0; k < stoi(input); k++)
	//{
	//	cout<< "(" << points[k][0] << ", " << points[k][1] << ")\n";
	//}

	//create an initial closestPair
	int closestPair [2] = {0, 1};

	//calculate the closestPair for each pair of Points
	//if the new distance is smaller than the old replace it

	auto t1 = std::chrono::high_resolution_clock::now();
	for(int i = 0; i < stoi(input); i++){
		for(int k = i+1; k < stoi(input); k++) {
			//old x and y values

			//dont need to calculate this everytime we go through loop might consider
			//taking it out
			long oldX = points[closestPair[0]][0]-points[closestPair[1]][0];
			long oldY = points[closestPair[0]][1]-points[closestPair[1]][1];
			long oldDistance = sqrt(pow(oldX,2)+pow(oldY,2));
			//new x and y values
			long newX = points[i][0]-points[k][0];
			long newY = points[i][1]-points[k][1];
			long newDistance = sqrt(pow(newX,2)+pow(newY,2));
			if(newDistance < oldDistance){
				closestPair[0] = i;
				closestPair[1] = k;
			}
		}
	}
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> fp_ms = t2 - t1;

	//prints the closestPair
	//cout << points[closestPair[0]][0] << "," << points[closestPair[0]][1] << " and " << points[closestPair[1]][0] << "," << points [closestPair[1]][1] << endl;
	cout << "it took: " << fp_ms.count() << "ms to complete" << endl;
}

void closest_grid(const char* input){
	long realLength = stol(input);
	cout << realLength << endl;

	int points [realLength][2];

	int closestPairGrid [3][3];
	for(long k = 0; k < stoi(input); k++)
	{
		points[k][0] = rand() %1000000;
		points[k][1] = rand() %1000000;
	}

	//calculate the max for x and y, done in linear time

	int xmax = points[0][0];
	int ymax = points[0][1];

	for(int i = 1; i < stoi(input); i++)
	{
		if(points[i][0]>xmax){
			xmax = points[i][0];

		}
		if(points[i][1]>ymax){
			ymax = points[i][1];
		}
	}

	xmax = (xmax/1000+1)*1000;
	ymax = (ymax/1000+1)*1000;


	//cout << "\n xmax is: " << xmax << ", " << "ymax is: " << ymax << "\n" << endl;

	int numPoints = stoi(input);
	int dimension = int(sqrt(numPoints));

	// vector<pair<int,int>> grid [dimension][dimension];
	vector<pair<int,int>>** grid = new vector<pair<int,int>>*[dimension];
	for(int i = 0; i < dimension; i++)
	{
		grid[i] = new vector<pair<int,int>>[dimension];
	}



	for(int i = 0; i < numPoints - 2; i++)
	{
		int xIndex = int(long(points[i][0])/long(xmax/dimension));
		if(xIndex == dimension){
			xIndex = xIndex -1;
		}
		int yIndex = int(long(points[i][1])/long(ymax/dimension));
		grid[xIndex][yIndex].push_back(make_pair(points[i][0],points[i][1]));


	}
	int xIndex = int(long(points[numPoints-2][0])/long(xmax/dimension));
	if(xIndex == dimension){
		xIndex = xIndex -1;
	}
	int yIndex = int(long(points[numPoints-2][1])/long(ymax/dimension));
	grid[xIndex][yIndex].push_back(make_pair(points[numPoints -2][0],points[numPoints -2][1]));

	int vecLength = grid[xIndex][yIndex].size();

	closestPairGrid[0][0] = xIndex;
	closestPairGrid[0][1] = yIndex;
	closestPairGrid[0][2] = vecLength;

	xIndex = int(long(points[numPoints-1][0])/long(xmax/dimension));
	if(xIndex == dimension){
		xIndex = xIndex -1;
	}
	yIndex = int(long(points[numPoints-1][1])/long(ymax/dimension));
	grid[xIndex][yIndex].push_back(make_pair(points[numPoints-1][0],points[numPoints-1][1]));

	vecLength = grid[xIndex][yIndex].size();


	closestPairGrid[1][0] = xIndex;
	closestPairGrid[1][1] = yIndex;
	closestPairGrid[1][2] = vecLength;

	long tempX = points[numPoints-1][0] - points[numPoints-2][0];
	long tempY = points[numPoints-1][1] - points[numPoints-2][1];
	long tempDist = sqrt(pow(tempX,2)+pow(tempY,2));


	closestPairGrid[2][0] = tempDist;
	// for(int i = 0; i < dimension; i ++){
	// 	for(int j = 0; j < dimension; j++){
	// 		cout << "index (i,j): " << i << ", " << j << endl;
	// 		for(int k = 0; k < grid[i][j].size(); k++){
	// 			cout << grid[i][j].at(k).first <<", " << grid[i][j].at(k).second << endl;
	// 		}
	// 		cout << "\n";
	// 	}
	// }
	auto t1 = std::chrono::high_resolution_clock::now();
	for( int i = 0; i < dimension; i++) {
		for(int j = 0; j< dimension; j++){

			if(grid[i][j].size() != 0){

				for( int k = 0; k < grid[i][j].size(); k ++){
					for(int n = k+1; n < grid[i][j].size(); n++){
						long newX = grid[i][j].at(k).first - grid[i][j].at(n).first;
						long newY = grid[i][j].at(k).second - grid[i][j].at(n).second;
						long newDist = sqrt(pow(newX,2)+pow(newY,2));
						if(newDist < closestPairGrid[2][0])
						{
							closestPairGrid[0][0] = i;
							closestPairGrid[0][1] = j;
							closestPairGrid[1][0] = i;
							closestPairGrid[1][1] = j;
							closestPairGrid[0][2] = k;
							closestPairGrid[1][2] = n;
							closestPairGrid[2][0] = newDist;
						}
					}
				}
			}
		}
	}

	for(int i = 0; i < dimension-1; i++)
	{
		for(int j = 0; j < dimension-1; j++)
		{

			if(grid[i][j].size() != 0)
			{
				for(int k = 0; k < grid[i][j].size(); k++)
				{
					for(int q = 0; q < grid[i+1][j].size(); q++)
					{
						long newX = grid[i][j].at(k).first - grid[i+1][j].at(q).first;
						long newY = grid[i][j].at(k).second - grid[i+1][j].at(q).second;
						long newDist = sqrt(pow(newX,2)+pow(newY,2));
						if(newDist < closestPairGrid[2][0])
						{
							closestPairGrid[0][0] = i;
							closestPairGrid[0][1] = j;
							closestPairGrid[1][0] = i+1;
							closestPairGrid[1][1] = j;
							closestPairGrid[0][2] = k;
							closestPairGrid[1][2] = q;
							closestPairGrid[2][0] = newDist;
						}
					}

					for(int q = 0; q < grid[i+1][j+1].size(); q++)
					{
						long newX = grid[i][j].at(k).first - grid[i+1][j+1].at(q).first;
						long newY = grid[i][j].at(k).second - grid[i+1][j+1].at(q).second;
						long newDist = sqrt(pow(newX,2)+pow(newY,2));
						if(newDist < closestPairGrid[2][0])
						{
							closestPairGrid[0][0] = i;
							closestPairGrid[0][1] = j;
							closestPairGrid[1][0] = i+1;
							closestPairGrid[1][1] = j+1;
							closestPairGrid[0][2] = k;
							closestPairGrid[1][2] = q;
							closestPairGrid[2][0] = newDist;
						}
					}

					for(int q = 0; q < grid[i][j+1].size(); q++)
					{
						long newX = grid[i][j].at(k).first - grid[i][j+1].at(q).first;
						long newY = grid[i][j].at(k).second - grid[i][j+1].at(q).second;
						long newDist = sqrt(pow(newX,2)+pow(newY,2));
						if(newDist < closestPairGrid[2][0])
						{

							closestPairGrid[0][0] = i;
							closestPairGrid[0][1] = j;
							closestPairGrid[1][0] = i;
							closestPairGrid[1][1] = j+1;
							closestPairGrid[0][2] = k;
							closestPairGrid[1][2] = q;
							closestPairGrid[2][0] = newDist;
						}
					}

				}
			}
		}
	}

	for(int i = 0; i < dimension-1; i++){
		int j = dimension-1;
		if(grid[i][j].size() != 0){
			for(int k = 0; k < grid[i][j].size(); k ++){
				for(int n = 0; n < grid[i+1][j].size(); n++){

					long newX = grid[i][j].at(k).first - grid[i+1][j].at(n).first;
					long newY = grid[i][j].at(k).second - grid[i+1][j].at(n).second;
					long newDist = sqrt(pow(newX,2)+pow(newY,2));
					if(newDist < closestPairGrid[2][0]){
						closestPairGrid[0][0] = i;
						closestPairGrid[0][1] = j;
						closestPairGrid[1][0] = i+1;
						closestPairGrid[1][1] = j;
						closestPairGrid[0][2] = k;
						closestPairGrid[1][2] = n;
						closestPairGrid[2][0] = newDist;
					}
				}
			}
		}


		if(grid[j][i].size() != 0){
			for( int k = 0; k < grid[j][i].size(); k ++){
				for(int n = 0; n < grid[j][i+1].size(); n++){
					long newX = grid[j][i].at(k).first - grid[j][i+1].at(n).first;
					long newY = grid[j][i].at(k).second - grid[j][i+1].at(n).second;
					long newDist = sqrt(pow(newX,2)+pow(newY,2));
					if(newDist < closestPairGrid[2][0])
					{
						closestPairGrid[0][0] = j;
						closestPairGrid[0][1] = i;
						closestPairGrid[1][0] = j;
						closestPairGrid[1][1] = i+1;
						closestPairGrid[0][2] = k;
						closestPairGrid[1][2] = n;
						closestPairGrid[2][0] = newDist;
					}
				}
			}
		}
	}
	auto t2 = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
	cout << "it took: " << fp_ms.count() << "ms to complete" << endl;
}


int main(int argc, const char * argv[]){
	srand(time(NULL)); //seeds the random number generation based on current time, is very random

	//closest_naive(argv[1]);
	closest_grid(argv[1]);

	return 0;

}
