//
//  RenderWindow.h
//  Renderer
//
//  Created by Raphael Kargon on 6/6/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#ifndef __Renderer__RenderWindow__
#define __Renderer__RenderWindow__

#include <QtWidgets/QFormLayout>
#include <QtWidgets/QWidget>
#include "RenderArea.h"

class RenderWindow : public QWidget
{
    //Q_OBJECT
    
public:
    RenderWindow();
private:
    RenderArea *renderArea;
};

#endif /* defined(__Renderer__RenderWindow__) */
