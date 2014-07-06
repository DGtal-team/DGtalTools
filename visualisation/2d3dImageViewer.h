/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/
/**
 * @file 2d3dImageViewer.h
 * @author Bertrand Kerautret (\c kerautre@loria.fr )
 * LORIA (CNRS, UMR 7503), University of Nancy, France
 *
 * @date 2014/07/04
 *
 * An example file named qglViewer.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#ifndef Q_MOC_RUN 
   #include "DGtal/io/viewers/Viewer3D.h"
   #include <DGtal/images/ImageContainerBySTLVector.h>
   #include "DGtal/helpers/StdDefs.h"
#endif 

#include <QMainWindow>

namespace Ui {
class MainWindow;

}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
explicit MainWindow(DGtal::Viewer3D<> *viewer,  DGtal::ImageContainerBySTLVector < DGtal::Z3i::Domain, unsigned char > *myImage3D,
                      QWidget *parent = 0, Qt::WindowFlags flags=0);
    ~MainWindow();
  void setImageProjX(const QPixmap &aPixMap);
  void setImageProjY(const QPixmap &aPixMap);
  void setImageProjZ(const QPixmap &aPixMap);


void updateSliceImageX(unsigned int sliceNumber, bool init);
void updateSliceImageY(unsigned int sliceNumber, bool init);
void updateSliceImageZ(unsigned int sliceNumber, bool init);

public slots:
  void updateSliceImageX();
  void updateSliceImageY();
  void updateSliceImageZ();


private:
    Ui::MainWindow *ui;
    DGtal::Viewer3D<> *myViewer;
    DGtal::ImageContainerBySTLVector < DGtal::Z3i::Domain, unsigned char > *myImage3D;
};

#endif // MAINWINDOW_H
  
