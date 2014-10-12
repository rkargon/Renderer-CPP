//
//  RenderWindow.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/6/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "RenderWindow.h"

RenderWindow::RenderWindow(){
    renderArea = new RenderArea(this);
    QLayout *layout = new QFormLayout();
    layout->addWidget(renderArea);
    this->setLayout(layout);
}