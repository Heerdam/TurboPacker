
#include <TP.hpp>

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/matrix_decompose.hpp>

#include <format>

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

    Detail::BoxList<float> postpacs;
    postpacs.ShuffleBoxes = true;
    postpacs.List.emplace_back(glm::vec3{28.f, 17.4f, 10.f}, 100);
    postpacs.List.emplace_back(glm::vec3{35.5f, 24.f, 12.5f}, 100);
    postpacs.List.emplace_back(glm::vec3{38.f, 35.f, 16.9f}, 100);
    postpacs.List.emplace_back(glm::vec3{53.5f, 28.5f, 16.5f}, 100);
    postpacs.List.emplace_back(glm::vec3{39.0f, 13.0f, 11.0f}, 100);
    postpacs.List.emplace_back(glm::vec3{55.5f, 37.0f, 6.0f}, 100);

    Config<float, CostFunction::CF_Krass> conf;
    //conf.MultiThreading = false;
    conf.NumThreads = 8;
    conf.UseRandomSeed_ = false;
    conf.Bounds = {80., 120.}; 
    conf.Height = 120.;
    conf.BoxType = Detail::BoxGenerationType::RANDOM;
    conf.CubeRandomBoxes = false;
    conf.LookAheadSize = 50;
    conf.EmptryTries = 0;
    conf.MaxEmptryTries = 0;
    conf.AllowOverlap = false;
    conf.MinBoxVolume = 15 * 15 * 15;
    conf.MaxBoxVolume = 30 * 30 * 30;
    conf.BoxList = std::move(postpacs);

    auto pr = solve(conf);
    pr.wait();

    const auto& b = pr.data();

    std::cout << std::endl << std::format("Done: [{}][{}%][{}s][Missed: {}]", pr.getBoxCount(), pr.getPackDensity() * 100., pr.getTime(), pr.getMissedCount()) << std::endl;

    //-------------------------------------

    Camera3D camera;
    camera.position = Vector3{ 0., 0., 250. };
    camera.target = Vector3{ 50.0f, 50.0f, 0.0f };  
    camera.up = Vector3{ 0.0f, 0.0f, 1.0f };
    camera.fovy = 45.0f;                       
    camera.projection = CAMERA_PERSPECTIVE;  

    float camDist = 250.;  // how far away from the target the camera is (radius)
    float rotAngle = 45; // the rotation angle around the target  (around Y)
    float tiltAngle = 5; // the tilt tangle of the camera (up/down)

    float rotSpeed = 0.25f; // to scale the mouse input
    float moveSpeed = 20.f; // to scale the linear input

    Vector2 cursorPos = GetMousePosition();

    //Shader shader = LoadShader(0, TextFormat("resources/shaders/glsl%i/grayscale.fs", 330));

    while (!WindowShouldClose()) {

         if (IsMouseButtonDown(1)) {
            Vector2 newPos = GetMousePosition();

            // update the angles from the delta
            rotAngle += (newPos.x - cursorPos.x) * rotSpeed;
            tiltAngle += (newPos.y - cursorPos.y) * rotSpeed;

            // clamp the tilt so we don't get gymbal lock
            if (tiltAngle > 89)
                tiltAngle = 89;
            if (tiltAngle < 1)
                tiltAngle = 1;
        }
        // always update the position so we don't get jumps
        cursorPos = GetMousePosition();

        // vector in rotation space to move
        Vector3 moveVec = { 0,0,0 };

        if (IsKeyDown(KEY_W))
            moveVec.z = -moveSpeed * GetFrameTime();
        if (IsKeyDown(KEY_S))
            moveVec.z = moveSpeed * GetFrameTime();

        if (IsKeyDown(KEY_A))
            moveVec.x = -moveSpeed * GetFrameTime();
        if (IsKeyDown(KEY_D))
            moveVec.x = moveSpeed * GetFrameTime();
    
        // update zoom
        camDist += GetMouseWheelMove();
        if (camDist < 1)
            camDist = 1;

        // vector we are going to transform to get the camera offset from the target point
        Vector3 camPos = { 0, 0, camDist};

        Matrix tiltMat = MatrixRotateX(tiltAngle * GetFrameTime()); // a matrix for the tilt rotation
        Matrix rotMat = MatrixRotateY(rotAngle * GetFrameTime()); // a matrix for the plane rotation
        Matrix mat = MatrixMultiply(tiltMat, rotMat); // the combined transformation matrix for the camera position

        camPos = Vector3Transform(camPos, mat); // transform the camera position into a vector in world space
        moveVec = Vector3Transform(moveVec, rotMat); // transform the movement vector into world space, but ignore the tilt so it is in plane

        camera.target = Vector3Add(camera.target, moveVec); // move the target to the moved position

        camera.position = Vector3Add(camera.target, camPos); // offset the camera position by the vector from the target positio

        //------------------------

        BeginDrawing();
        ClearBackground(RAYWHITE);

        BeginMode3D(camera);
        //DrawGrid(25, 1.0f);

        const auto ext = conf.Bounds * 0.5f;
        DrawCubeWiresV(Vector3{ext.y, ext.x, conf.Height * 0.5f}, Vector3{conf.Bounds.y, conf.Bounds.x, (float)conf.Height}, BLUE);

        for (const auto& tr : b) {
            const auto trans = glm::vec<3, float>(tr[3]);
            DrawCubeV(Vector3{ trans.x, trans.y, trans.z }, Vector3{ 
                glm::length(glm::vec3(tr[0])), 
                glm::length(glm::vec3(tr[1])), 
                glm::length(glm::vec3(tr[2])) 
                }, RED);
            DrawCubeWiresV(Vector3{ trans.x, trans.y, trans.z }, Vector3{ 
                glm::length(glm::vec3(tr[0])), 
                glm::length(glm::vec3(tr[1])), 
                glm::length(glm::vec3(tr[2])) 
                }, BLUE);
        }

        EndMode3D();
        EndDrawing();

    }

    CloseWindow();

    return 0;
}