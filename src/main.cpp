
#include <TP.hpp>

#include <glm/gtc/type_ptr.hpp>

#include <raylib.h>
#include <raymath.h>

void taskgraphtest() {

    using namespace TP;
    using namespace TP::Detail;

    std::mutex m;
    int32_t i = 0;

    const auto f = [&]() {
        std::lock_guard<std::mutex> lock(m);
        std::cout << i++ << std::endl;
    };

    TaskGraph tg (4);

    for(int32_t i = 0; i < 10; ++i)
        tg.dispatch(f);

    tg.wait();

}

int main() {

    using namespace TP;

    InitWindow(1000, 1000, "TurboPacker");
    SetTargetFPS(60); 

    Config<float, CostFunction::CF_Basic> conf;
    conf.Bounds = {40., 50. };
    conf.Height = 60.;
    conf.EmptryTries = 1;
    conf.MinBoxVolume = 20 * 20 * 20;
    conf.MaxBoxVolume = 30 * 30 * 30;

    const auto pr = solve(conf);

    while(!pr.isDone()){
        std::cout << pr.getBoxCount() << " [" << pr.getPackDensity() << "]"  << "\r";
    }


    const auto& b = pr.data();
    
    std::cout << std::endl << "done [" << b.size() << "]" << std::endl;


    //-------------------------------------

    Camera camera;
    camera.position = Vector3{ 5., 5., 5. };
    camera.target = Vector3{ 0.0f, 0.0f, 0.0f };  
    camera.up = Vector3{ 0.0f, 1.0f, 0.0f };
    camera.fovy = 45.0f;                       
    camera.projection = CAMERA_PERSPECTIVE;  


    while (!WindowShouldClose()) {

        BeginDrawing();
        ClearBackground(RAYWHITE);

        BeginMode3D(camera);
        DrawGrid(25, 1.0f);

        for (const auto& tr : b) {
            const Matrix& mat = *reinterpret_cast<const Matrix*>(glm::value_ptr(tr));
            DrawCubeWiresV(Vector3{ mat.m12, mat.m13, mat.m14 }, Vector3{ 
                glm::length(glm::vec3(tr[0])), 
                glm::length(glm::vec3(tr[1])),  
                glm::length(glm::vec3(tr[2])), 
                }, RED);
        }

        EndMode3D();
        EndDrawing();

    }

    CloseWindow();

    return 0;
}