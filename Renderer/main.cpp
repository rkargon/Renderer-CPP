//
//  main.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/5/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include <iostream>
#include <fstream>

#include <QtWidgets/QApplication>
#include <QtWidgets/QgraphicsScene>
#include <QtWidgets/QgraphicsView>
#include <QtWidgets/QMainWindow>

#include "RenderWindow.h"

using namespace std;

int main(int argc, char **argv)
{
    //srand((unsigned int)time(0));
    QApplication app(argc, argv);
    RenderWindow win;
    win.resize(400, 400);
    win.show();
    
    return app.exec();
}

