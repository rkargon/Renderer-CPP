//
//  View.cpp
//  Renderer
//
//  Created by Raphael Kargon on 6/6/14.
//  Copyright (c) 2014 Raphael Kargon. All rights reserved.
//

#include "view.h"

#include "glm/gtx/io.hpp"

#include <QApplication>
#include <QCommandLineParser>
#include <QFileInfo>

struct view_options {
  render_options render_opts;
  std::string scene_file;
};

view_options process_cli_args() {
  QCommandLineParser parser;
  parser.setApplicationDescription("Test helper");
  parser.addHelpOption();

  parser.addOption({"threads", "Number of render threads. A value of 0 means "
                               "to auto-detect thread num.",
                    "threads", "0"});
  parser.addOption(
      {"samples", "Number of path tracing samples", "samples", "10"});
  parser.addOption({"file", "STL or OBJ file to load", "file"});
  parser.addOption({"depth", "Max number of ray bounces", "depth", "10"});
  parser.addOption({"tilesize", "Tile size", "tilesize", "32"});
  parser.addOption({"antialiasing_per_side", "Antialiasing samples per side",
                    "antialiasing_per_side", "1"});
  parser.process(*QApplication::instance());

  return {{parser.value("threads").toUInt(), parser.value("samples").toUInt(),
           parser.value("depth").toUInt(), parser.value("tilesize").toUInt(),
           parser.value("antialiasing_per_side").toUInt()},
          parser.value("file").toStdString()};
}

View::View(QWidget *parent) : QWidget(parent), imgrasters(width(), height()) {
  setBackgroundRole(QPalette::Base);
  setAutoFillBackground(true);
  setFocusPolicy(Qt::StrongFocus);
  setFocus(Qt::ActiveWindowFocusReason);
  setLayout(new QHBoxLayout);

  auto opts = process_cli_args();

  renderimg = QImage(reinterpret_cast<uchar *>(imgrasters.colbuffer.get()),
                     width(), height(), QImage::Format_ARGB32);
  this->init_scene(opts.scene_file);
  manager =
      std::make_unique<thread_manager>(&imgrasters, &sc, opts.render_opts);

  // status label
  statuslbl.reset(new QLabel("Raph Renderer 2018"));
  statuslbl->setParent(this);
  this->layout()->addWidget(statuslbl.get());
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
    for (auto &o_ptr : sc.objects) {
      o_ptr->smooth = !o_ptr->smooth;
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

  for (const auto &obj_ptr : sc.objects) {
    for (const edge &e : obj_ptr->edges) {
      p1 = sc.cam.project_vertex(obj_ptr->vertices[e[0]], w, h);
      p2 = sc.cam.project_vertex(obj_ptr->vertices[e[1]], w, h);
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

void View::init_scene(const std::string &scene_file) {
  std::string models_dir =
      "/Users/raphaelkargon/Documents/Programming/Computer Graphics/Models/";
  std::ifstream dragonfile(models_dir + "dragonsmall.stl");
  std::ifstream spherefile(models_dir + "sphere.stl");
  std::ifstream cubefile(models_dir + "cube.stl");

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
    sc.w = std::make_unique<world>();
    sc.materials.emplace_back(color(1));
    sc.bsdfs.emplace_back(new EmissionBSDF());
    sc.bsdfs.emplace_back(new DiffuseBSDF());

    mesh &sphere = sc.add_object(spherefile, "Sphere", false);
    sphere.mat_id = 0;
    sphere.bsdf = sc.bsdfs[0].get();
    // sphere.project_texture(TEX_PROJ_SPHERICAL);
    sphere.scale_centered(vertex(0.3, 0.3, 0.3));
    sphere.move(vertex(0, 0, 2));

    mesh &dragon = sc.add_object(dragonfile, "Dragon", false);
    dragon.mat_id = 0;
    dragon.bsdf = sc.bsdfs[1].get();
    // dragon.project_texture(TEX_PROJ_SPHERICAL);
    sc.update_tree();

    sc.de_obj = de_mandelbox(2.5, 1, 0.5, 1, 30);
  } else {

    ///////
    sc = scene(QFileInfo(QString::fromStdString(scene_file))
                   .absoluteFilePath()
                   .toStdString());
    sc.w = std::make_unique<world>();
    sc.lamps.emplace_back(100, 2, vertex(10, 10, 10), rgb_to_color(0xFFFFFF));
    sc.lamps.emplace_back(100, 2, vertex(10, 10, 10), rgb_to_color(0xFFFFFF));
    sc.lamps.emplace_back(100, 2, vertex(10, 10, 10), rgb_to_color(0xFFFFFF));
    sc.de_obj = de_mandelbox(2.5, 1, 0.5, 1, 30);
    sc.cam.center = vertex(0, .8, 3.6);
    sc.cam.focus = vertex(0);
    sc.cam.normal = vertex(0, 0, -1);
    sc.cam.vert = vertex(0, 1, 0);
    sc.cam.calc_image_vectors();

    for (const auto &obj_ptr : sc.objects) {
      std::cout << *obj_ptr << std::endl;
      //  std::cout << sc.materials[obj_ptr->mat_id].name << std::endl;
      std::cout << sc.materials[obj_ptr->mat_id].diff_col << std::endl;
    }
    sc.update_tree();
  }
}
