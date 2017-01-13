
#include <iostream>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/registration/icp.h>
#include <pcl/features/normal_3d.h>

#include <boost/thread/thread.hpp>
#include <pcl/common/common_headers.h>
#include <pcl/visualization/pcl_visualizer.h>

#include <pcl/filters/voxel_grid.h>

int main(int argc, char** argv)
{
	// set to 1 to apply previous transformation to subsequent data
	int rolling_transform = 1;
	double xT = 0, yT = 0, zT = 0;
	double R11T = 1, R12T = 0, R13T = 0;
	double R21T = 0, R22T = 1, R23T = 0;
	double R31T = 0, R32T = 0, R33T = 1;
	Eigen::Matrix4f transform_new;

	double occupancy_voxels = 0;
	float voxelGridResolution = 0.05;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			if (i == j) {
				transform_new(i, j) = 1;
			}
			else
			{
				transform_new(i, j) = 0;
			}
		}
	}




	// input file strings
	std::string fixedFilename, queryFilename;

	// uGPS data sizes
	int N_hits = 541;						// number of SICK laser scan points

	// declare point clouds
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_A_in(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_B_in(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_A(new pcl::PointCloud<pcl::PointNormal>);
	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_B(new pcl::PointCloud<pcl::PointNormal>);
	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_A_sub(new pcl::PointCloud<pcl::PointNormal>);
	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_B_sub(new pcl::PointCloud<pcl::PointNormal>);
	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_B_subT(new pcl::PointCloud<pcl::PointNormal>);

	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_A_sub_F(new pcl::PointCloud<pcl::PointNormal>);
	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_B_sub_F(new pcl::PointCloud<pcl::PointNormal>);
	pcl::PointCloud<pcl::PointNormal>::Ptr cloud_B_subT_F(new pcl::PointCloud<pcl::PointNormal>);


	pcl::PointCloud<pcl::Normal>::Ptr cloud_A_normals(new pcl::PointCloud<pcl::Normal>);
	pcl::PointCloud<pcl::Normal>::Ptr cloud_B_normals(new pcl::PointCloud<pcl::Normal>);

	// ask user for input cloud filenames
	std::cout << "Enter fixed cloud filename:" << std::endl;
	std::cin >> fixedFilename;
	std::cout << "Enter query cloud filename:" << std::endl;
	std::cin >> queryFilename;

	std::cout << "Enter number of laser hits per scan slice: ";
	std::cin >> N_hits;

	// ICP parameters
	int N_wnd_B;			// query data window size
	int N_wnd_A;  // fixed data window size (specified in terms of query data indices, will be converted later)

	// ask user for data window sizes
	std::cout << "Enter query data window size: ";
	std::cin >> N_wnd_B;
	std::cout << "Enter fixed data window size (wrt query window indices): ";
	std::cin >> N_wnd_A;

	double r_normals;
	int N_closest = 10;

	std::cout << "Enter number of closest points for normals: ";
	std::cin >> N_closest;
	std::cout << "Apply rolling transform (1 for yes)? ";
	std::cin >> rolling_transform;


	//std::cout << "Enter radius for computing normals: ";
	//std::cin >> r_normals;


	// load pcd files
	std::cout << "Loading files..." << std::endl;

	// fixed cloud
	pcl::io::loadPCDFile<pcl::PointXYZ>(fixedFilename, *cloud_A_in);
	int N_pts_A = cloud_A_in->points.size();						// number of points in cloud A
	int N_scans_A = N_pts_A / N_hits;							// number of SICK laser scan slices (rings) in cloud A
	std::cout << "		fixed cloud loaded: " << N_scans_A << " scans | " << N_pts_A << " points" << std::endl;

	// query cloud
	pcl::io::loadPCDFile<pcl::PointXYZ>(queryFilename, *cloud_B_in);
	int N_pts_B = cloud_B_in->points.size();						// number of points in cloud B
	int N_scans_B = N_pts_B / N_hits;							// number of SICK laser scan slices (rings) in cloud B
	std::cout << "		query cloud loaded: " << N_scans_B << " scans | " << N_pts_B << " points" << std::endl;


	//------------------------------------------------------------------------------------
	std::cout << "Pre-processing..." << std::endl;

	// fill cloud_A data		
	cloud_A->width = N_pts_A;
	cloud_A->height = 1;
	cloud_A->is_dense = false;
	cloud_A->points.resize(N_pts_A);

	// fill cloud_B data		
	cloud_B->width = N_pts_B;
	cloud_B->height = 1;
	cloud_B->is_dense = false;
	cloud_B->points.resize(N_pts_B);
	
	int nanFlagA2 = 0;
	int nanFlagB2 = 0;

	for (int j = 0; j < N_pts_A; j++)
	{
		cloud_A->points[j].x = cloud_A_in->points[j].x;
		cloud_A->points[j].y = cloud_A_in->points[j].y;
		cloud_A->points[j].z = cloud_A_in->points[j].z;
	}

	for (int j = 0; j < N_pts_B; j++)
	{
		cloud_B->points[j].x = cloud_B_in->points[j].x;
		cloud_B->points[j].y = cloud_B_in->points[j].y;
		cloud_B->points[j].z = cloud_B_in->points[j].z;
	}

	
	std::cout << "      Computing normals..." << std::endl;


	int nanFlagA = 1;
	int nanFlagB = 1;
	
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normalEstimator;
	pcl::search::KdTree<pcl::PointXYZ>::Ptr searchTree(new pcl::search::KdTree<pcl::PointXYZ>());
	normalEstimator.setSearchMethod(searchTree);


	while ((nanFlagA != 0) || (nanFlagB != 0)) {

		//normalEstimator.setRadiusSearch(r_normals);
		normalEstimator.setKSearch(N_closest);

		normalEstimator.setInputCloud(cloud_A_in);
		normalEstimator.compute(*cloud_A_normals);

		normalEstimator.setInputCloud(cloud_B_in);
		normalEstimator.compute(*cloud_B_normals);

		nanFlagA = 0;
		nanFlagB = 0;
		
		for (int j = 0; j < N_pts_A; j++)
		{
			cloud_A->points[j].normal_x = cloud_A_normals->points[j].normal_x;
			cloud_A->points[j].normal_y = cloud_A_normals->points[j].normal_y;
			cloud_A->points[j].normal_z = cloud_A_normals->points[j].normal_z;

			if (std::isnan(cloud_A->points[j].normal_x)) {
				nanFlagA = 1;
			}
		}

		for (int j = 0; j < N_pts_B; j++)
		{
			cloud_B->points[j].normal_x = cloud_B_normals->points[j].normal_x;
			cloud_B->points[j].normal_y = cloud_B_normals->points[j].normal_y;
			cloud_B->points[j].normal_z = cloud_B_normals->points[j].normal_z;

			if (std::isnan(cloud_B->points[j].normal_x)) {
				nanFlagB = 1;
			}
		}

		

		if ((nanFlagA == 1) || (nanFlagB == 1)) {

			N_closest++;

			std::cout << " Cloud normals contain NaN elements, increasing N_closest to " << N_closest << " -- nanFlagA = " << nanFlagA << " -- nanFlagB = " << nanFlagB << endl;
		}
	}

	/*
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);

	//... read, pass in or create a point cloud ...

	// Create the normal estimation class, and pass the input dataset to it
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> ne;
	ne.setInputCloud(cloud);

	// Create an empty kdtree representation, and pass it to the normal estimation object.
	// Its content will be filled inside the object, based on the given input dataset (as no other search surface is given).
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());
	ne.setSearchMethod(tree);

	// Output datasets
	pcl::PointCloud<pcl::Normal>::Ptr cloud_normals(new pcl::PointCloud<pcl::Normal>);

	// Use all neighbors in a sphere of radius 3cm
	ne.setRadiusSearch(0.03);

	// Compute the features
	ne.compute(*cloud_normals);

	*/



	std::cout << "		Computing centroids..." << std::endl;
	// find scan slice centroids
	pcl::PointCloud<pcl::PointXYZ>::Ptr centroids_A(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr centroids_B(new pcl::PointCloud<pcl::PointXYZ>);
	int k;

	centroids_A->width = N_scans_A;
	centroids_A->height = 1;
	centroids_A->is_dense = false;
	centroids_A->points.resize(N_scans_A);

	centroids_B->width = N_scans_B;
	centroids_B->height = 1;
	centroids_B->is_dense = false;
	centroids_B->points.resize(N_scans_B);

	std::cout << "		...fixed" << std::endl;
	// fixed data
	for (int i = 0; i < N_scans_A; i++)
	{
		centroids_A->points[i].x = 0;
		centroids_A->points[i].y = 0;
		centroids_A->points[i].z = 0;

		for (int j = 0; j < N_hits; j++)
		{
			k = i*N_hits + j;
			centroids_A->points[i].x += cloud_A->points[k].x;
			centroids_A->points[i].y += cloud_A->points[k].y;
			centroids_A->points[i].z += cloud_A->points[k].z;
		}
		centroids_A->points[i].x /= N_hits;
		centroids_A->points[i].y /= N_hits;
		centroids_A->points[i].z /= N_hits;
	}
	std::cout << "			...fixed centroids computed" << std::endl;

	std::cout << "		...query" << std::endl;
	// query data
	for (int i = 0; i < N_scans_B; i++)
	{
		centroids_B->points[i].x = 0;
		centroids_B->points[i].y = 0;
		centroids_B->points[i].z = 0;
		for (int j = 0; j < N_hits; j++)
		{
			k = i*N_hits + j;
			centroids_B->points[i].x += cloud_B->points[k].x;
			centroids_B->points[i].y += cloud_B->points[k].y;
			centroids_B->points[i].z += cloud_B->points[k].z;
		}
		centroids_B->points[i].x /= N_hits;
		centroids_B->points[i].y /= N_hits;
		centroids_B->points[i].z /= N_hits;
	}
	std::cout << "			...query centroids computed" << std::endl;

	std::cout << "Correcting normal directions" << endl;

	float cpx, cpy, cpz, nx, ny, nz, px, py, pz;
	int flipA = 0;
	int flipB = 0;

	for (int i = 0; i < N_scans_A; i++)
	{
		
		cpx = centroids_A->points[i].x;
		cpy = centroids_A->points[i].y;
		cpz = centroids_A->points[i].z;

		for (int j = 0; j < N_hits; j++) {

			k = i*N_hits + j;

			px = cloud_A->points[k].x;
			py = cloud_A->points[k].y;
			pz = cloud_A->points[k].z;

			nx = cloud_A->points[k].normal_x;
			ny = cloud_A->points[k].normal_y;
			nz = cloud_A->points[k].normal_z;

			if (((cpx - px)*nx + (cpy - py)*ny + (cpz - pz)*nz) < 0) {
				cloud_A->points[k].normal_x = -nx;
				cloud_A->points[k].normal_y = -ny;
				cloud_A->points[k].normal_z = -nz;
				flipA++;
			}
		}	
	}

	for (int i = 0; i < N_scans_B; i++)
	{

		cpx = centroids_B->points[i].x;
		cpy = centroids_B->points[i].y;
		cpz = centroids_B->points[i].z;

		for (int j = 0; j < N_hits; j++) {

			k = i*N_hits + j;

			px = cloud_B->points[k].x;
			py = cloud_B->points[k].y;
			pz = cloud_B->points[k].z;

			nx = cloud_B->points[k].normal_x;
			ny = cloud_B->points[k].normal_y;
			nz = cloud_B->points[k].normal_z;

			if (((cpx - px)*nx + (cpy - py)*ny + (cpz - pz)*nz) < 0) {
				cloud_B->points[k].normal_x = -nx;
				cloud_B->points[k].normal_y = -ny;
				cloud_B->points[k].normal_z = -nz;
				flipB++;
			}
		}
	}

	cout << "flipA = " << flipA << endl;
	cout << "flipB = " << flipB << endl;


	// output new pcd files with normals
	pcl::io::savePCDFileASCII("output_fixed_cloud.pcd", *cloud_A);
	pcl::io::savePCDFileASCII("output_query_cloud.pcd", *cloud_B);
	pcl::io::savePCDFileASCII("output_fixed_cloud_centroids.pcd", *centroids_A);
	pcl::io::savePCDFileASCII("output_query_cloud_centroids.pcd", *centroids_B);

	
	std::cout << "		Finding closest centroids" << std::endl;

	// find closest scan based on centroids (ie. list of closest fixed scan for each query scan)
	double centroid_distance, min_distance;
	int j_min;
	std::vector<int> closest_idx(N_scans_B);

	for (int i = 0; i < N_scans_B; i++)
	{
		min_distance = 100000.0;
		for (int j = 0; j < N_scans_A; j++)
		{
			centroid_distance = sqrt(pow(centroids_A->points[j].x - centroids_B->points[i].x, 2) + pow(centroids_A->points[j].y - centroids_B->points[i].y, 2) + pow(centroids_A->points[j].z - centroids_B->points[i].z, 2));
			if (centroid_distance < min_distance)
			{
				j_min = j;
				min_distance = centroid_distance;
			}
		}
		closest_idx[i] = j_min;
		//if (i < 50)
		//{
		//	std::cout << "			" << i << " - " << j_min << std::endl;
		//}
	}
	std::cout << "			...closest centroids computed" << std::endl;

	/*
	// opening viewer
	boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("3D Viewer"));
	viewer->setBackgroundColor(0, 0, 0);
	viewer->addPointCloud<pcl::PointXYZ>(cloud_A_in, "reference cloud");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 3, "reference cloud");
	viewer->addPointCloudNormals<pcl::PointXYZ, pcl::Normal>(cloud_A_in, cloud_A_normals, 10, 0.05, "normals");
	viewer->addCoordinateSystem(1.0);
	viewer->initCameraParameters();
	
	while (!viewer->wasStopped())
	{
		viewer->spinOnce(100);
		boost::this_thread::sleep(boost::posix_time::microseconds(100000));
	}
	*/


	//=====================================================================================
	// MOVING WINDOW ICP SOLUTION
	
	// create ICP structure
	pcl::IterativeClosestPointWithNormals<pcl::PointNormal, pcl::PointNormal> icp;
	pcl::PointCloud<pcl::PointNormal> Final;
	Eigen::Matrix4f transform;

	

	
	std::cout << "Process windowed ICP solutions" << std::endl;

	// declare variables
	int i_min_B, i_max_B;
	int k_min_B, k_max_B;
	int k_min_A, k_max_A;
	int N_A, N_B;
	int i_min_AB, i_max_AB, i_min_A, i_max_A;
	float mean_fitness;

	

	// setup output file
	std::string fnameOut = "output_temp.txt";
	FILE* fileOut = fopen(fnameOut.c_str(), "w");


	// loop through each slice in query scan
	for (int i = 0; i < N_scans_B; i++)
	{
		// set indices for query window

		i_min_B = i - N_wnd_B;
		if (i_min_B < 0)
			i_min_B = 0;
		i_max_B = i + N_wnd_B;
		if (i_max_B >= N_scans_B)
			i_max_B = N_scans_B - 1;



		k_min_B = i_min_B*N_hits;
		k_max_B = i_max_B*N_hits - 1;

		N_B = k_max_B - k_min_B + 1;

		// fill cloud_B_sub data		
		cloud_B_sub->width = N_B;
		cloud_B_sub->height = 1;
		cloud_B_sub->is_dense = false;
		cloud_B_sub->points.resize(N_B);

		for (int j = 0; j < N_B; j++)
		{
			cloud_B_sub->points[j].x = cloud_B->points[k_min_B + j].x;
			cloud_B_sub->points[j].y = cloud_B->points[k_min_B + j].y;
			cloud_B_sub->points[j].z = cloud_B->points[k_min_B + j].z;
			cloud_B_sub->points[j].normal_x = cloud_B->points[k_min_B + j].normal_x;
			cloud_B_sub->points[j].normal_y = cloud_B->points[k_min_B + j].normal_y;
			cloud_B_sub->points[j].normal_z = cloud_B->points[k_min_B + j].normal_z;
		}



		// fill cloud_A_sub data
		i_min_AB = i - N_wnd_A;
		i_max_AB = i + N_wnd_A;

		if (i_min_AB < 0)
			i_min_AB = 0;

		if (i_max_AB >= N_scans_B)
			i_max_AB = N_scans_B - 1;

		i_min_A = closest_idx[i_min_AB];
		i_max_A = closest_idx[i_max_AB];

		int i_temp;
		if (i_min_A > i_max_A)
		{
			i_temp = i_min_A;
			i_min_A = i_max_A;
			i_max_A = i_temp;
		}

		k_min_A = i_min_A*N_hits;
		k_max_A = i_max_A*N_hits - 1;

		N_A = k_max_A - k_min_A + 1;

		cloud_A_sub->width = N_A;
		cloud_A_sub->height = 1;
		cloud_A_sub->is_dense = false;
		cloud_A_sub->points.resize(N_A);

		for (int j = 0; j < N_A; j++)
		{
			cloud_A_sub->points[j].x = cloud_A->points[k_min_A + j].x;
			cloud_A_sub->points[j].y = cloud_A->points[k_min_A + j].y;
			cloud_A_sub->points[j].z = cloud_A->points[k_min_A + j].z;
			cloud_A_sub->points[j].normal_x = cloud_A->points[k_min_A + j].normal_x;
			cloud_A_sub->points[j].normal_y = cloud_A->points[k_min_A + j].normal_y;
			cloud_A_sub->points[j].normal_z = cloud_A->points[k_min_A + j].normal_z;
		}

		//std::cout << N_A << ", " << N_B << " | " << i_min_A << ", " << i_max_A << " | " << i_min_B << ", " << i_max_B << " | " << k_min_A << ", " << k_max_A << " | " << k_min_B << ", " << k_max_B << std::endl;

		std::cout << i_min_A << "-" << i_max_A << " | " << N_A << " --- " << i_min_B << "-" << i_max_B << " | " << N_B << std::endl;

		if (rolling_transform == 1) {
			pcl::transformPointCloud(*cloud_B_sub, *cloud_B_subT, transform_new);
		}

		
		if (occupancy_voxels == 1) {
			// down-sample clouds to set resolution
			pcl::VoxelGrid<pcl::PointNormal> voxelGridDownSample;
			voxelGridDownSample.setInputCloud(cloud_A_sub);
			voxelGridDownSample.setLeafSize(voxelGridResolution, voxelGridResolution, voxelGridResolution);
			voxelGridDownSample.filter(*cloud_A_sub_F);

			voxelGridDownSample.setInputCloud(cloud_B_sub);
			voxelGridDownSample.setLeafSize(voxelGridResolution, voxelGridResolution, voxelGridResolution);
			voxelGridDownSample.filter(*cloud_B_sub_F);

			voxelGridDownSample.setInputCloud(cloud_B_subT);
			voxelGridDownSample.setLeafSize(voxelGridResolution, voxelGridResolution, voxelGridResolution);
			voxelGridDownSample.filter(*cloud_B_subT_F);
			
			//======================================================================================
			// run ICP

			icp.setInputTarget(cloud_A_sub_F);
			if (rolling_transform == 1) {
				icp.setInputSource(cloud_B_subT_F);
			}
			else {
				icp.setInputSource(cloud_B_sub_F);
			}
		}
		else {
			//======================================================================================
			// run ICP

			icp.setInputTarget(cloud_A_sub);
			if (rolling_transform == 1) {
				icp.setInputSource(cloud_B_subT);
			}
			else {
				icp.setInputSource(cloud_B_sub);
			}
		}
		
		

		icp.align(Final);

		// get final transformation
		transform = icp.getFinalTransformation();
		transform_new = transform*transform_new;

		// write to file... [px py pz R11 R12 R13 R21 R22 R23 R31 R32 R33 fitness]
		for (int j = 0; j < 3; j++)
			fprintf(fileOut, "%f ", transform(j, 3));

		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
				fprintf(fileOut, "%f ", transform(j, k));
		}
		fprintf(fileOut, "%f ", icp.getFitnessScore());
		fprintf(fileOut, "\n");

		if (i == 0)
			mean_fitness = icp.getFitnessScore();
		else
			mean_fitness = (i*mean_fitness + icp.getFitnessScore()) / (i + 1);

		std::cout << "Scan " << i << " |  " << icp.hasConverged() << " | " << icp.getFitnessScore() << " | " << mean_fitness << std::endl;

	}

	fclose(fileOut);

	cout << "Press ENTER to exit..." << endl;
	cin.sync();
	cin.get();

	return 0;
}