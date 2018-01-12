#ifndef BEMTOOLTESTS_MISC_VIEW_HPP
#define BEMTOOLTESTS_MISC_VIEW_HPP

#include <nanogui/nanogui.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <iostream>
#include <fstream>

typedef std::complex<double> K;

#include <htool/view.hpp>
#include <bemtool/tools.hpp>

namespace bemtool{

template<int dim>
std::vector<htool::N4> translate_elts(const Dof<BasisFct<P1,dim>>& dof){
    int nb_elt = NbElt(dof);
    std::vector<htool::N4> Elts(nb_elt);
    for (int i=0;i<nb_elt;i++){
        for (int j=0;j<dim+1;j++){
            Elts[i][j]=dof[i][j];
        }
    }
    return Elts;
}

template<int dim>
std::vector<htool::R3> translate_normal(const Mesh<dim>& mesh){
    int nb_elt = NbElt(mesh);
    std::vector<htool::R3> Normals(nb_elt);
    const std::vector<R3>& normals = NormalTo(mesh);
    for (int i=0;i<nb_elt;i++){
        for (int j=0;j<3;j++){
            Normals[i][j]=normals[i][j];
        }
    }
    return Normals;
}

void attach_ui(htool::Scene& s) {
    htool::statics& gv = htool::Scene::gv;

    nanogui::Window *window = new nanogui::Window(htool::Scene::gv.screen, "BemTool");
    window->setPosition(Eigen::Vector2i(650, 200));
    window->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Vertical,nanogui::Alignment::Middle, 10, 10));
    nanogui::Widget* tools = new nanogui::Widget(window);
    tools->setLayout(new nanogui::BoxLayout(nanogui::Orientation::Vertical, nanogui::Alignment::Middle, 0, 6));

    nanogui::ComboBox* combobox_eq  = new nanogui::ComboBox(window, { "LA", "HE", "YU"});
    nanogui::ComboBox* combobox_op  = new nanogui::ComboBox(window, { "SL", "DL", "TDL", "HS"});


    nanogui::Button* b = new nanogui::Button(tools, "Set operator");
    b->setCallback([&,combobox_eq,combobox_op] {
        if (gv.active_project == NULL)
        std::cerr << "No active project" << std::endl;
        else{
            std::string str = nanogui::file_dialog({{"msh", "Mesh file"}}, false);
            glfwFocusWindow(gv.glwindow);
            if (str!=""){
                if (combobox_op->selectedIndex()==0 && combobox_eq->selectedIndex()==0){
                    double kappa = 1;
                    // Mesh
                    Geometry* node = new Geometry(str);
                    Mesh<2>* mesh= new Mesh<2>();
                    mesh->Load(*node,0);
                    Orienting(*mesh);
                    int nb_elt = NbElt(*mesh);

                    // Dof
                    Dof<P1_2D> dof(*mesh);
                    int nb_dof = NbDof(dof);
                    std::vector<htool::R3> x(nb_dof);
                    for (int i=0;i<nb_dof;i++){
                        x[i][0]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][0];
                        x[i][1]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][1];
                        x[i][2]=dof(((dof.ToElt(i))[0])[0])[((dof.ToElt(i))[0])[1]][2];
                    }

                    // Generator
                    SubBIO_Generator<LA_SL_3D_P1xP1,P1_2D>* generator = new SubBIO_Generator<LA_SL_3D_P1xP1,P1_2D>(dof,kappa);



                    // GLMesh
                    std::vector<int> NbPts(nb_elt,3);
                    htool::GLMesh m(x,translate_elts(dof),NbPts,translate_normal(*mesh));
                    std::vector<int> tab(x.size());
                    std::iota(tab.begin(),tab.end(),int(0));
                    m.set_tab(tab);

                    // Set project
                    s.set_mesh(m);
                    gv.active_project->set_matrix(generator);
                    gv.active_project->set_ctrs(x);
                    gv.active_project->set_rays(std::vector<double> (x.size(),0));
                }

            }
        }
    });

    // b = new nanogui::Button(tools, "Load matrix");
    // b->setCallback([&] {
    //     if (gv.active_project == NULL)
    //     std::cerr << "No active project" << std::endl;
    //     else{
    //         std::string strmat = nanogui::file_dialog(
    //             {{"bin", "Matrix binary file"}}, false);
    //             std::cout << "Loading matrix file " << strmat << " ..." << std::endl;
    //             htool::Matrix<K> *A = new htool::Matrix<K>;
    //             A->bytes_to_matrix(strmat);
    //             gv.active_project->set_matrix(A);
    //     }
    // });

    gv.screen->performLayout();

}

} // namespace
#endif
