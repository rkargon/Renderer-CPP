//
//  main.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/5/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "ui/mainwindow.h"

#include <QApplication>

#include <iostream>

int main(int argc, char **argv) {
  QApplication app(argc, argv);
  MainWindow w;
  w.resize(400, 400);
  w.show();

  return app.exec();
}
