
#include <TP.hpp>

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/matrix_decompose.hpp>

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
    conf.Bounds = {80., 100. };
    conf.Height = 30.;
    conf.EmptryTries = 100;
    conf.MinBoxVolume = 15 * 15 * 15;
    conf.MaxBoxVolume = 30 * 30 * 30;

    const auto pr = solve(conf);

    while(!pr.isDone()){
        //std::osyncstream(std::cout) << pr.getBoxCount() << " [" << pr.getPackDensity() << "]"  << "\r";
    }


    const auto& b = pr.data();
    //for (const auto& tr : b)
        //std::cout << glm::to_string(tr) << std::endl;
    
    std::cout << std::endl << "done [" << b.size() << "]" << std::endl;


    //-------------------------------------

    Camera camera;
    camera.position = Vector3{ 150., 150., 150. };
    camera.target = Vector3{ 0.0f, 0.0f, 0.0f };  
    camera.up = Vector3{ 0.0f, 0.0f, 1.0f };
    camera.fovy = 45.0f;                       
    camera.projection = CAMERA_PERSPECTIVE;  


    while (!WindowShouldClose()) {

        BeginDrawing();
        ClearBackground(RAYWHITE);

        BeginMode3D(camera);
        //DrawGrid(25, 1.0f);

        const auto ext = conf.Bounds * 0.5f;
        DrawCubeWiresV(Vector3{ext.y, ext.x, conf.Height * 0.5f}, Vector3{conf.Bounds.y, conf.Bounds.x, (float)conf.Height}, BLUE);

        for (const auto& tr : b) {
            const auto trans = glm::vec<3, float>(tr[3]);
            DrawCubeWiresV(Vector3{ trans.x, trans.y, trans.z }, Vector3{ 
                glm::length(glm::vec3(tr[0])), 
                glm::length(glm::vec3(tr[1])), 
                glm::length(glm::vec3(tr[2])) 
                }, RED);
        }

        EndMode3D();
        EndDrawing();

    }

    CloseWindow();

    return 0;
}