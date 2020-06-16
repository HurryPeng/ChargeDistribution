#ifndef _HURRYPENG_DEBUGUTIL_HPP
#define _HURRYPENG_DEBUGUTIL_HPP

#include "Vector3D.hpp"
#include "ElectricField.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <list>
#include <vector>

void exportElectricField(std::string filename, const HurryPeng::ElectricField & electricField, int scale = 50)
{
    std::ofstream ofs(filename);
    ofs << "{\n";
    for (int i = -scale; i < scale; i++) for (int j = -scale; j < scale; j++)
    {
        HurryPeng::Vector3D temp = electricField.get({double(i) / 1000, double(j) / 1000, 0});
        ofs << "    {{" << double(i) / 1000 << ", " << double(j) / 1000 << "}, {" << temp.x << ", " << temp.y << "}}";
        if (i != scale - 1 || j != scale - 1) ofs << ',';
        ofs << '\n';
    }
    ofs << "}\n";
    ofs.close();
}

template <typename Container>
void exportPointSet(std::string filename, const Container & points)
{
    std::ofstream ofs(filename);
    ofs << "{\n";
    bool firstLine = true;
    for (auto & point : points)
    {
        if (!firstLine) ofs << ", \n";
        else firstLine = false;
        ofs << "    " << point;
        auto it = points.end();
        it--;
    }
    ofs << "\n}\n";
    ofs.close();
}

template <typename Container>
void exportPointSetForPy(std::string filename, const Container & points)
{
    std::ofstream ofs(filename);
    ofs << "[\n";
    bool firstLine = true;
    for (auto & point : points)
    {
        if (!firstLine) ofs << ", \n";
        else firstLine = false;
        ofs << "    " << "[" << point.x << ", " << point.y << ", " << point.z << "]";
        auto it = points.end();
        it--;
    }
    ofs << "\n]\n";
    ofs.close();
}

template <typename Container>
std::list<HurryPeng::Vector3D> chargesToPointSet(const Container & charges)
{
    std::list<HurryPeng::Vector3D> pointSet;
    for (const HurryPeng::FreeCharge & charge : charges)
        pointSet.push_back(charge.coord);
    return pointSet;
}

std::vector<long double> paramStatU(std::vector<std::vector<std::optional<long double>>> fields)
// fields[uInt][vInt] = intensity
{
    int precision = fields.size();
    std::vector<long double> stat(precision);
    for (int uInt = 0; uInt < precision; uInt++)
    {
        long double sum = 0;
        int count = 0;
        for (const std::optional<long double> & intensity : fields[uInt])
        {
            if (intensity)
            {
                sum += *intensity;
                count++;
            }
        }
        stat[uInt] = sum / count;
    }
    return stat;
}

#endif
