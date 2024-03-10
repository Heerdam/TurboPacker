
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

    InitWindow(1000, 1000, "Test");
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

    std::cout << "done" << std::endl;

    const auto& b = pr.data();

    //-------------------------------------

    Camera camera;
    camera.position = Vector3{ 100., 100., 100. };
    camera.target = Vector3{ 0.0f, 0.0f, 0.0f };  
    camera.up = Vector3{ 0.0f, 0.0f, 1.0f };
    camera.fovy = 45.0f;                       
    camera.projection = CAMERA_PERSPECTIVE;  


    while (!WindowShouldClose()) {

        BeginMode3D(camera);

        ClearBackground(RAYWHITE);

        DrawGrid(10, 1.0f);

        for (const auto& tr : b) {
            const Matrix& mat = *reinterpret_cast<const Matrix*>(glm::value_ptr(tr));
            DrawCubeV(Vector3{ mat.m12, mat.m13, mat.m14 }, Vector3{ 
                glm::length(glm::vec3(tr[0])), 
                glm::length(glm::vec3(tr[1])),  
                glm::length(glm::vec3(tr[2])), 
                }, RED);
        }

        EndMode3D();

    }

    CloseWindow();

    return 0;
}