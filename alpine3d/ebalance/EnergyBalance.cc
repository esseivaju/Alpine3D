/***********************************************************************************/
/*  Copyright 2009-2015 WSL Institute for Snow and Avalanche Research    SLF-DAVOS      */
/***********************************************************************************/
/* This file is part of Alpine3D.
    Alpine3D is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Alpine3D is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with Alpine3D.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <alpine3d/ebalance/EnergyBalance.h>
#include <alpine3d/MPIControl.h>
#include <alpine3d/OMPControl.h>

using namespace mio;
using namespace std;

// FELIX: added pv_pts, pv_points
EnergyBalance::EnergyBalance(const unsigned int& i_nbworkers, const mio::Config& cfg_in, const mio::DEMObject &dem_in)
              : snowpack(NULL), terrain_radiation(NULL), radfields(i_nbworkers), dem(dem_in), vecMeteo(),
                albedo(dem_in, 0.), direct_unshaded_horizontal(), direct(), diffuse(), reflected(),
                timer(), dimx(dem_in.getNx()), dimy(dem_in.getNy()), nbworkers(i_nbworkers), pv_points(), PVP(nullptr), cfg(cfg_in)
{

	MPIControl& instance = MPIControl::instance();

	size_t startx = 0, nx = dimx;
	instance.getArraySliceParams(dimx, startx, nx);

	#pragma omp parallel for schedule(static)
	for (size_t ii=0; ii<nbworkers; ii++) {
		size_t thread_startx, thread_nx;
		OMPControl::getArraySliceParams(nx, nbworkers, ii, thread_startx, thread_nx);
		const size_t offset = startx + thread_startx;

		#pragma omp critical(ebWorkers_status)
		std::cout << "[i] EnergyBalance worker " << ii << " on process " << instance.rank() << " will start at offset " << offset << " with nx " << thread_nx << "\n";
		radfields[ii] = RadiationField(dem_in, offset, thread_nx);
	}

	if (instance.master())
		std::cout << "[i] EnergyBalance initialized a total of " << instance.size() << " process(es) with " << nbworkers << " worker(s) each\n";

	if (cfg.keyExists("PVPFILE", "EBalance"))
	{
		//load PVP data
		readPVP();
		//performs some validation on loaded data
		std::vector<Coords> co_vec;
		for (size_t ii = 0; ii < pv_points.size(); ii++)
		{
			Coords point;
			point.setXY(pv_points[ii][0], pv_points[ii][1], pv_points[ii][2]);
			co_vec.push_back(point);
		}
		bool master = instance.master();
		if (!dem.gridify(co_vec, true))
		{ //keep invalid points
			if (master)
				cerr << "[E] Some PVP are invalid or outside the DEM:\n";
			for (size_t ii = 0; ii < co_vec.size(); ii++)
				if (!co_vec[ii].indexIsValid() && master)
					std::cout << "[E] Point " << ii << "\t" << co_vec[ii].toString(Coords::CARTESIAN) << "\n";
			throw InvalidArgumentException("Invalid PVP, please check in the logs", AT);
		}
		else if (master)
			std::cout << "[i] Using " << pv_points.size() << " PVP\n";

		if(cfg.get("Terrain_Radiation", "EBalance")) PVP= new SolarPanel(cfg, dem, pv_points, &radfields[0]);
	}


	// Every MPI process will have its own copy of terrain_radiation object with full DEM
	const bool enable_terrain_radiation = cfg.get("Terrain_Radiation", "EBalance");
	if (enable_terrain_radiation) {
		terrain_radiation = TerrainRadiationFactory::getAlgorithm(cfg, dem, nbworkers, PVP);
		const std::string algo = terrain_radiation->algo;
		if (instance.master())
			std::cout << "[i] Using terrain radiation with model: " << algo << "\n";
	}
}

EnergyBalance::~EnergyBalance() {
	Destroy( );
}

EnergyBalance& EnergyBalance::operator=(const EnergyBalance& source) {
	if (this != &source) {
		snowpack = source.snowpack;
		terrain_radiation = source.terrain_radiation;
		radfields = source.radfields;
		dem = source.dem;
		vecMeteo = source.vecMeteo;
		albedo = source.albedo;
		direct = source.direct;
		diffuse = source.diffuse;
		reflected = source.reflected;
		direct_unshaded_horizontal = source.reflected;
		timer = source.timer;
		dimx = source.dimx;
		dimy = source.dimy;
		nbworkers = source.nbworkers;
	}
	return *this;
}

std::string EnergyBalance::getGridsRequirements() const
{
	return "TOP_ALB";
}

void EnergyBalance::Destroy()
{
	if (terrain_radiation) {
		delete terrain_radiation;
		terrain_radiation = NULL;
	}
	if (PVP) {
		delete PVP;
		PVP = nullptr;
	}
}

void EnergyBalance::setSnowPack(SnowpackInterface& mysnowpack)
{
	snowpack = &mysnowpack;
}

void EnergyBalance::setAlbedo(const mio::Grid2DObject& in_albedo)
{
	albedo = in_albedo;

	direct_unshaded_horizontal.resize(0, 0); //FELIX
	direct.resize(0, 0); //resetting these grids that are not valid anymore
	diffuse.resize(0, 0);
	reflected.resize(0, 0);
}

void EnergyBalance::setStations(const std::vector<mio::MeteoData>& in_vecMeteo)
{
	vecMeteo = in_vecMeteo;

	direct_unshaded_horizontal.resize(0, 0); //FELIX
	direct.resize(0, 0); //resetting these grids that are not valid anymore
	diffuse.resize(0, 0);
	reflected.resize(0, 0);
}

void EnergyBalance::setMeteo(const mio::Grid2DObject& in_ilwr,
                             const mio::Grid2DObject& in_ta, const mio::Grid2DObject& in_rh, const mio::Grid2DObject& in_p, const mio::Date timestamp)
{
	timer.restart();
	direct.resize(dimx, dimy);
	diffuse.resize(dimx, dimy);
	direct_unshaded_horizontal.resize(dimx, dimy);

	#pragma omp parallel for schedule(dynamic)
	for (size_t ii=0; ii<nbworkers; ii++) {
		radfields[ii].setStations(vecMeteo, albedo); //calculate the parameters at the radiation stations
		size_t startx, nx;
		radfields[ii].getBandOffsets(startx, nx);
		radfields[ii].setMeteo(mio::Grid2DObject(in_ta, startx, 0, nx, dimy),
		                       mio::Grid2DObject(in_rh, startx, 0, nx, dimy),
		                       mio::Grid2DObject(in_p, startx, 0, nx, dimy),
		                       mio::Grid2DObject(albedo, startx, 0, nx, dimy));

		mio::Array2D<double> band_direct, band_diffuse, band_direct_unshaded_horizontal;
		radfields[ii].getRadiation(band_direct, band_diffuse, band_direct_unshaded_horizontal);
		direct.fill(band_direct, startx, 0, nx, dimy);
		diffuse.fill(band_diffuse, startx, 0, nx, dimy);
		direct_unshaded_horizontal.fill(band_direct_unshaded_horizontal, startx, 0, nx, dimy);
	}
	MPIControl::instance().allreduce_sum(direct);
	MPIControl::instance().allreduce_sum(diffuse);
	MPIControl::instance().allreduce_sum(direct_unshaded_horizontal);
	double solarAzimuth, solarElevation;
	radfields[0].getPositionSun(solarAzimuth, solarElevation);

	if (cfg.keyExists("PVPFILE", "EBalance"))
	{
		PVP->setGridRadiation(albedo, direct, diffuse, direct_unshaded_horizontal);
	}

	mio::Array2D<double> view_factor(dimx, dimy, IOUtils::nodata);
	if (terrain_radiation) {
		// note: parallelization has to take place inside the TerrainRadiationAlgorithm implementations
		terrain_radiation->setMeteo(albedo.grid2D, in_ta.grid2D, in_rh.grid2D, in_ilwr.grid2D);
		terrain_radiation->getRadiation(direct, diffuse, reflected, direct_unshaded_horizontal,view_factor,
                                    solarAzimuth, solarElevation);
	}

	if (MPIControl::instance().master())
		cout << "[i] Ebalance simulation done for " << timestamp.toString(Date::ISO) << "\n";

	if (snowpack) {
		double solarAzimuth, solarElevation;
		radfields[0].getPositionSun(solarAzimuth, solarElevation); //we need it only for handing over to snowpack

		mio::Array2D<double> ilwr = in_ilwr.grid2D;
		mio::Array2D<double> global = direct+diffuse; //otherwise the compiler does not match the types

		if (!reflected.empty()) global += reflected;

		timer.stop();
		try {
			snowpack->setRadiationComponents(global, ilwr, diffuse, view_factor, reflected,
                                       in_ilwr.grid2D,  solarElevation, timestamp); //this triggers Snowpack calculation
		} catch(std::exception& e) {
			std::cout << "[E] Exception in snowpack->setRadiationComponents()\n";
			cout << e.what() << endl;
			std::abort(); //force core dump
		}
	}
	timer.stop();
}

void EnergyBalance::setPVP(const mio::Date timestamp){

	if (cfg.keyExists("PVPFILE", "EBalance")) PVP->setPVP(timestamp);
}

void EnergyBalance::writeSumPVP(const unsigned int max_steps){

	if (cfg.keyExists("PVPFILE", "EBalance")) PVP->writeSumPVP(max_steps);
}

void EnergyBalance::readPVP()
{
	const std::string filename = cfg.get("PVPFILE", "EBalance");
	if (!FileUtils::fileExists(filename))
	{
		throw NotFoundException(filename, AT);
	}

	smet::SMETReader myreader(filename);
	std::vector<double> vec_data;
	myreader.read(vec_data);
	const size_t nr_fields = myreader.get_nr_of_fields();
	const int epsg = myreader.get_header_intvalue("epsg");
	const double smet_nodata = myreader.get_header_doublevalue("nodata");

	if (myreader.location_in_data(smet::WGS84) == true)
	{
		size_t lat_fd = IOUtils::unodata, lon_fd = IOUtils::unodata, wid_fd = IOUtils::unodata, hig_fd = IOUtils::unodata;
		size_t alt_fd = IOUtils::unodata, inc_fd = IOUtils::unodata, az_fd = IOUtils::unodata;
		for (size_t ii = 0; ii < nr_fields; ii++)
		{
			const std::string tmp(myreader.get_field_name(ii));
			if (tmp == "latitude")
				lat_fd = ii;
			if (tmp == "longitude")
				lon_fd = ii;
			if (tmp == "altitude")
				alt_fd = ii;
			if (tmp == "inclination")
				inc_fd = ii;
			if (tmp == "azimuth")
				az_fd = ii;
			if (tmp == "height")
				hig_fd = ii;
			if (tmp == "width")
				wid_fd = ii;
		}
		for (size_t ii = 0; ii < vec_data.size(); ii += nr_fields)
		{
			std::vector<double> point;
			point = {vec_data[ii + lat_fd], vec_data[ii + lon_fd], vec_data[ii + alt_fd], vec_data[ii + inc_fd], vec_data[ii + az_fd], vec_data[ii + hig_fd], vec_data[ii + wid_fd]};
			pv_points.push_back(point);
		}
	}
	else if (myreader.location_in_data(smet::EPSG) == true)
	{
		if (epsg == (int)floor(smet_nodata + 0.1))
			throw InvalidFormatException("In file \"" + filename + "\", missing EPSG code in header!", AT);

		size_t east_fd = IOUtils::unodata, north_fd = IOUtils::unodata, alt_fd = IOUtils::unodata;
		size_t inc_fd = IOUtils::unodata, az_fd = IOUtils::unodata, wid_fd = IOUtils::unodata, hig_fd = IOUtils::unodata;
		for (size_t ii = 0; ii < nr_fields; ii++)
		{
			const std::string tmp(myreader.get_field_name(ii));
			if (tmp == "easting")
				east_fd = ii;
			if (tmp == "northing")
				north_fd = ii;
			if (tmp == "altitude")
				alt_fd = ii;
			if (tmp == "inclination")
				inc_fd = ii;
			if (tmp == "azimuth")
				az_fd = ii;
			if (tmp == "height")
				hig_fd = ii;
			if (tmp == "width")
				wid_fd = ii;
		}
		if ((east_fd == IOUtils::unodata) || (north_fd == IOUtils::unodata) || (alt_fd == IOUtils::unodata))
			throw InvalidFormatException("File \"" + filename + "\" does not contain all data fields necessary for EPSG coordinates", AT);

		for (size_t ii = 0; ii < vec_data.size(); ii += nr_fields)
		{
			Coords coord_temp;
			coord_temp.setEPSG(epsg);
			coord_temp.setXY(vec_data[ii + east_fd], vec_data[ii + north_fd], vec_data[ii + alt_fd]);

			std::vector<double> point;
			//point={coord_temp.getLat(), coord_temp.getLon(), vec_data[ii+alt_fd], vec_data[ii+inc_fd], vec_data[ii+az_fd]};
			point = {coord_temp.getEasting(), coord_temp.getNorthing(), vec_data[ii + alt_fd], vec_data[ii + inc_fd], vec_data[ii + az_fd], vec_data[ii + hig_fd], vec_data[ii + wid_fd]};
			pv_points.push_back(point);
		}
	}
	else
	{
		throw InvalidFormatException("File \"" + filename + "\" does not contain expected location information in DATA section!", AT);
	}
}

double EnergyBalance::getTiming() const
{
	return timer.getElapsed();
}
