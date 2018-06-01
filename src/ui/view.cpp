//
//  View.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/6/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "view.h"

View::View(QWidget *parent) : QWidget(parent), imgrasters(width(), height()) {
  setBackgroundRole(QPalette::Base);
  setAutoFillBackground(true);
  setFocusPolicy(Qt::StrongFocus);
  setFocus(Qt::ActiveWindowFocusReason);
  setLayout(new QHBoxLayout);

  renderimg = QImage(reinterpret_cast<uchar *>(imgrasters.colbuffer.get()),
                     width(), height(), QImage::Format_ARGB32);
  this->init_scene();
  manager = new thread_manager(&imgrasters, &sc, 1, 0);

  // status label
  statuslbl = new QLabel("Raph Renderer 2018");
  statuslbl->setParent(this);
  this->layout()->addWidget(statuslbl);
  statuslbl->setSizePolicy(QSizePolicy::MinimumExpanding,
                           QSizePolicy::MinimumExpanding);
  statuslbl->setAlignment(Qt::AlignTop);
  updateText();
  statuslbl->show();
}

QSize View::minimumSizeHint() const { return QSize(400, 400); }

QSize View::sizeHint() const { return QSize(10000, 10000); }

void View::paintEvent(QPaintEvent *) {
  QPainter painter(this);
  painter.drawImage(QPoint(0, 0), renderimg);
}

void View::updateText() {
  switch (rendermode) {
  case 0:
    statuslbl->setText("1. Wireframe");
    break;
  case 1:
    statuslbl->setText("2. ZBuffer Draw");
    break;
  case 2:
    statuslbl->setText("3. SSAO");
    break;
  case 3:
    statuslbl->setText("4. Normal map");
    break;
  case 4:
    statuslbl->setText("5. Raytracing");
    break;
  case 5:
    statuslbl->setText("6. Path Tracing");
    break;
  case 6:
    statuslbl->setText("6. Ambient Occlusion");
    break;
  case 7:
    statuslbl->setText("7. Fractal Rendering");
    break;
  default:
    statuslbl->setText("Raph Renderer 2015");
  }
}

void View::updateImage() {
  manager->stop();
  switch (rendermode) {
  case 0:
    drawWireFrame();
    break;
  case 1:
    zbuffer_draw(imgrasters, sc);
    break;
  case 2:
    SSAO(imgrasters, sc);
    break;
  case 3:
    paint_normal_map(imgrasters, sc);
    break;
  case 4:
    manager->set_render_method(ray_trace_pixel);
    manager->start();
    break;
  case 5:
    manager->set_render_method(path_trace_pixel);
    manager->start();
    break;
  case 6:
    manager->set_render_method(amb_occ_pixel);
    manager->start();
    break;
  case 7:
    manager->set_render_method(ray_march_pixel);
    manager->start();
    break;
  }
}

/* INPUT LISTENERS */

void View::keyPressEvent(QKeyEvent *event) {
  switch (event->key()) {
  case Qt::Key_5:
    sc.cam.ortho = !sc.cam.ortho;
    updateImage();
    break;
  case Qt::Key_S:
    // TODO: add per-object smoothing and, in general, selection
    for (mesh &o : sc.objects) {
      o.smooth = !o.smooth;
    }
    updateText();
    updateImage();
    break;
  case Qt::Key_Z:
    if (event->modifiers() & Qt::ShiftModifier) {
      --rendermode;
    } else {
      ++rendermode;
    }
    rendermode = (rendermode + 8) % 8;
    updateText();
    updateImage();
    break;
  case Qt::Key_Space:
    sc.cam.focus = vertex();
    sc.cam.center_focus();
    updateImage();
    break;
  default:
    QWidget::keyPressEvent(event);
  }
  repaint();
}

void View::mouseMoveEvent(QMouseEvent *event) {
  QPoint pos = event->pos();
  QPoint delta = pos - prevpos; // prevpos initialized in mousePressEvent

  if (event->modifiers() & Qt::ShiftModifier) {
    sc.cam.shift_focus(-delta.x() / 100.0, delta.y() / 100.0);
  } else {
    sc.cam.rotate_local_x(delta.y() / 100.0);
    sc.cam.rotate_local_y(delta.x() / 100.0);
  }
  sc.cam.center_focus();
  prevpos = pos;
  updateImage();
  repaint();
}

void View::mousePressEvent(QMouseEvent *event) {
  prevpos = event->pos();
  if (event->button() == Qt::MiddleButton) {
    sc.cam.set_global_rotation(0, 0, 0);
    sc.cam.center_focus();
    updateImage();
    repaint();
  }
}

void View::wheelEvent(QWheelEvent *event) {
  double zoomfactor = pow(1.005, -event->delta());
  sc.cam.zoom((double)zoomfactor);
  updateImage();
  repaint();
}

void View::resizeEvent(QResizeEvent *event) {
  // thread manager will segfault if image buffer changes out from under its
  // feet.
  manager->stop();
  int w = event->size().width();
  int h = event->size().height();
  imgrasters.resize(w, h);
  renderimg = QImage(reinterpret_cast<uchar *>(imgrasters.colbuffer.get()), w,
                     h, QImage::Format_RGB32);
  updateImage();
}

/* DRAWING FUNCTIONS */

void View::drawWireFrame() {
  renderimg.fill(0xffffff);
  QPainter painter(&renderimg);
  int w = width(), h = height();
  point_2d<double> p1, p2;
  painter.setPen(QColor(0, 0, 0));

  for (const mesh &obj : sc.objects) {
    for (const edge &e : obj.edges) {
      p1 = sc.cam.project_vertex(obj.vertices[e[0]], w, h);
      p2 = sc.cam.project_vertex(obj.vertices[e[1]], w, h);
      if (isnan(p1.x) || isnan(p2.x) || isnan(p1.y) || isnan(p2.y)) {
        continue;
      }
      painter.drawLine(p1.x, p1.y, p2.x, p2.y);
    }
  }

  bool draw_kdtree = false;
  if (draw_kdtree) {
    painter.setPen(QColor(255, 0, 0));
    for (const auto &e : sc.kdt.wireframe()) {
      p1 = sc.cam.project_vertex(e.first, w, h);
      p2 = sc.cam.project_vertex(e.second, w, h);
      if (isnan(p1.x) || isnan(p2.x) || isnan(p1.y) || isnan(p2.y)) {
        continue;
      }
      painter.drawLine(p1.x, p1.y, p2.x, p2.y);
    }
  }
}

QColor colorToQColor(const color &c) {
  auto c1 = glm::clamp(c * 255.0, color(0.0), color(255.0));
  return QColor((int)c1.r, (int)c1.g, (int)c1.b);
}

void View::init_scene() {
  std::string models_dir =
      "/Users/raphaelkargon/Documents/Programming/Computer Graphics/Models/";
  std::ifstream dragonfile(models_dir + "dragonsmall.stl");
  std::ifstream spherefile(models_dir + "sphere.stl");
  std::ifstream cubefile(models_dir + "cube.stl");

#undef STL_SCENE
#ifdef STL_SCENE
  bool stl_scene = false;
  if (stl_scene) {
    // scene setup
    sc.cam.center = vertex(0, -5, 0);
    sc.cam.normal = vertex(0, 1, 0);
    sc.cam.focus = vertex(0, 0, 0);
    sc.cam.dof_focus_distance = 4;
    sc.cam.aperture_size = 0.0;
    sc.cam.center_focus();
    sc.cam.calc_image_vectors();
    sc.lamps.emplace_back(25, 2, vertex(0, 5, 5), rgb_to_color(0xFF66CC));
    sc.w = std::make_unique<sky>();

    mesh &sphere = sc.add_object(spherefile, "Sphere", false);
    sphere.mat = new material(color(1));
    sphere.bsdf = new EmissionBSDF();
    // sphere.project_texture(TEX_PROJ_SPHERICAL);
    sphere.scale_centered(vertex(0.3, 0.3, 0.3));
    sphere.move(vertex(0, 0, 1));

    mesh &dragon = sc.add_object(dragonfile, "Dragon", false);
    dragon.mat = new material(color(1));
    dragon.bsdf = new DiffuseBSDF();
    // dragon.project_texture(TEX_PROJ_SPHERICAL);
    sc.update_tree();

    sc.de_obj = de_mandelbox(2.5, 1, 0.5, 1, 30);
  }
#endif

  ///////
  sc = scene(models_dir + "cornell-box/CornellBox-Sphere.obj");
  sc.w = std::make_unique<sky>();
  sc.lamps.emplace_back(100, 2, vertex(0, 10, 10), rgb_to_color(0xFF66CC));
  sc.lamps.emplace_back(100, 2, vertex(10, 0, 10), rgb_to_color(0x66FFCC));
  sc.lamps.emplace_back(100, 2, vertex(10, 10, 00), rgb_to_color(0xCC66FF));
  sc.de_obj = de_mandelbox(2.5, 1, 0.5, 1, 30);
  sc.cam.center = vertex(0, .8, 3.6);
  sc.cam.focus = vertex(0);
  sc.cam.normal = vertex(0, 0, -1);
  sc.cam.vert = vertex(0, 1, 0);
  sc.cam.calc_image_vectors();

  for (const auto &obj : sc.objects) {
    std::cout << obj << std::endl;
  }
  sc.update_tree();
}
