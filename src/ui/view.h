//
//  View.h
//  Renderer
//
//  Created by Raphael Kargon on 6/6/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef _RENDERER_UI_VIEW_H
#define _RENDERER_UI_VIEW_H

#include "camera.h"
#include "distance_estimation.h"
#include "kdtree.h"
#include "mesh.h"
#include "raster.h"
#include "rendering.h"
#include "renderthreads.h"

#include <QtCore/QTimer>
#include <QtGui/QColor>
#include <QtGui/QKeyEvent>
#include <QtGui/QMouseEvent>
#include <QtGui/QPainter>
#include <QtGui/QPalette>
#include <QtGui/QWheelEvent>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLayout>
#include <QtWidgets/QWidget>

#include <algorithm>
#include <fstream>
#include <thread>

class View : public QWidget {
  Q_OBJECT

public:
  scene sc;

  unsigned int rendermode = 0;
  raster imgrasters;
  QImage renderimg;
  thread_manager *manager;

  QPoint prevpos; // for mouse movement tracking
  QLabel *statuslbl;

  View(QWidget *parent = 0);
  ~View() = default;
  QSize minimumSizeHint() const;
  QSize sizeHint() const;

  void updateText();
  void updateImage();
  void drawWireFrame();

protected:
  void paintEvent(QPaintEvent *event);
  void keyPressEvent(QKeyEvent *event);
  void mouseMoveEvent(QMouseEvent *event);
  void mousePressEvent(QMouseEvent *event);
  void wheelEvent(QWheelEvent *event);
  void resizeEvent(QResizeEvent *event);
};

QColor colorToQColor(const color &c);

#endif // _RENDERER_UI_VIEW_H
