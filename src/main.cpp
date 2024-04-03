
#include <TP.hpp>
#include <Util.hpp>

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <glm/gtx/matrix_decompose.hpp>

#include <format>
#include <charconv>
#include <future>

#include <raylib.h>
#include <raymath.h>

#include <imgui.h>
#include <imgui_impl_raylib.h>

#include <font_regular.h>
#include <font_bold.h>


int main() {

    using namespace TP;

    InitWindow(1920, 1080, "TurboPacker");

    const auto scale = GetWindowScaleDPI();

    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard; 
    ImFont* regular = io.Fonts->AddFontFromMemoryTTF(Montserrat_Regular_ttf, sizeof(Montserrat_Regular_ttf), 14.f * std::max(scale.x, scale.y));
    ImFont* bold = io.Fonts->AddFontFromMemoryTTF(Montserrat_Bold_ttf, sizeof(Montserrat_Bold_ttf), 16.f * std::max(scale.x, scale.y));
    ImGui::StyleColorsDark();
    
    ImGui_ImplRaylib_Init();
    Imgui_ImplRaylib_BuildFontAtlas();

    Detail::BoxList<float> postpacs;
    postpacs.ShuffleBoxes = true;
    postpacs.List.emplace_back(glm::vec3{28.f, 17.4f, 10.f}, 100);
    postpacs.List.emplace_back(glm::vec3{35.5f, 24.f, 12.5f}, 100);
    postpacs.List.emplace_back(glm::vec3{38.f, 35.f, 16.9f}, 100);
    postpacs.List.emplace_back(glm::vec3{53.5f, 28.5f, 16.5f}, 100);
    postpacs.List.emplace_back(glm::vec3{39.0f, 13.0f, 11.0f}, 100);
    postpacs.List.emplace_back(glm::vec3{55.5f, 37.0f, 6.0f}, 100);

    Config<float, CostFunction::CF_Krass> conf;
    conf.MultiThreading = true;
    conf.NumThreads = 8;
    conf.UseRandomSeed = true;
    conf.Seed = 12341234;
    conf.Bounds = {80., 120.}; 
    conf.Height = 120.;
    conf.BoxType = Detail::BoxGenerationType::VALIDATE;
    conf.CubeRandomBoxes = false;
    conf.LookAheadSize = 50;
    conf.EmptryTries = 0;
    conf.MaxEmptryTries = 0;
    conf.AllowOverlap = false;
    conf.MinBoxVolume = 15 * 15 * 15;
    conf.MaxBoxVolume = 30 * 30 * 30;
    conf.BoxList = std::move(postpacs);

    using Prmise = decltype(solve(conf));
    std::unique_ptr<Prmise> pr = nullptr;
    //-------------------------------------
    Util::Camera3d camera;
    camera.up = { 0., 1., 0. };
    camera.target = { conf.Bounds.x * 0.5f, 0., conf.Bounds.y * 0.5f };
    camera.camDist = 250.f;
    camera.tiltAngle = -65.f;
    //-------------------------------------

    std::vector<glm::mat<4, 4, float>> b;
    int32_t cc = 0;

    const auto col = [&](float _scalar) -> Color {
        const float RedSclr = std::clamp<float>((1.0f - _scalar)/0.5f,0.f,1.f);
        const float GreenSclr = std::clamp<float>((_scalar/0.5f),0.f,1.f);
        const uint8_t R = (uint8_t)(255 * RedSclr);
        const uint8_t G = (uint8_t)(255 * GreenSclr);
        const uint8_t B = 0;
        return Color(R, G, B, 255);
    };

    

    const char* str_modes[] = { "Random", "List", "Validate" };
    int cur_mode = (int)conf.BoxType;
    uint64_t seedval = conf.Seed;
    std::unique_ptr<Util::Detail::EvalRes> eval_res = nullptr;

    //-------------------------------------

    while (!WindowShouldClose()) {

        ImGui_ImplRaylib_ProcessEvents();
       
        //------------------------

        ImGui_ImplRaylib_NewFrame();
        ImGui::NewFrame();
        ImGui::SetNextWindowSize(ImVec2(475, 700), ImGuiCond_FirstUseEver);
        ImGui::Begin("PackerWidget", (bool*)nullptr, ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize);

        if(pr) {
            ImGui::PushFont(bold);
            ImGui::Text(std::format("{}s", pr->getTime()).data());
            ImGui::PopFont();
            ImGui::Dummy({0, 10});
            ImGui::Text("Progress");
            ImGui::ProgressBar(pr->getPackDensity());
            if(pr->isDone()){
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.f, 1.f, 0.f, 1.f));
                ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.0f, 0.0f, 0.0f, 1.0f));
                if(ImGui::Button("Pack", ImVec2(100, 30))) {
                    pr = std::make_unique<Prmise>(solve(conf));
                }
                ImGui::PopStyleColor(2);
            } else {
                ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(1.f, 0.f, 0.0f, 1.f));
                if(ImGui::Button("Stop", ImVec2(100, 30))) {
                    pr->stop();
                    pr = nullptr;
                }
                ImGui::PopStyleColor(1);
            }
        } else {
            ImGui::PushFont(bold);
            ImGui::Text(std::format("{}s", 0.).data());
            ImGui::PopFont();
            ImGui::Dummy({0, 10});
            ImGui::Text("Progress");
            ImGui::ProgressBar(0.f);
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.0f, 0.0f, 0.0f, 1.0f));
            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.f, 1.f, 0.f, 1.f));
            if(ImGui::Button("Pack", ImVec2(100, 30))) {

                
                if(conf.BoxType == Detail::BoxGenerationType::VALIDATE && eval_res != nullptr) {
                    Config<float, CostFunction::CF_Krass> cc = conf;
                    cc.BoxType = Detail::BoxGenerationType::LIST;

                    std::sort(eval_res->decomp_.begin(), eval_res->decomp_.end(), [](const auto& _v1, const auto& _v2){
                        return _v1.first < _v2.first;
                    });

                    Detail::BoxList<float> list;
                    list.ShuffleBoxes = true;
                    for(const auto& [id, bx] : eval_res->decomp_){
                        list.List.push_back( { {bx.getSize().x, bx.getSize().y, conf.Height - 1.f}, 1});
                    }
                    cc.BoxList = list;

                    pr = std::make_unique<Prmise>(solve(cc));
                } else  pr = std::make_unique<Prmise>(solve(conf));
               
            }
            ImGui::PopStyleColor(2);
        }

        ImGui::Dummy({0, 10});
        ImGui::Separator();
        ImGui::Dummy({0, 4});
        ImGui::PushFont(bold);
        ImGui::Text("Config");
        ImGui::PopFont();
        ImGui::Dummy({0, 4});

        if(ImGui::BeginCombo("Mode", str_modes[cur_mode])) {
            for(int32_t i = 0; i < IM_ARRAYSIZE(str_modes); ++i){
                if(ImGui::Selectable(str_modes[i], cur_mode == i)){
                    cur_mode = i;
                    switch(cur_mode){
                        case 0: conf.BoxType = Detail::BoxGenerationType::RANDOM; break;
                        case 1: conf.BoxType = Detail::BoxGenerationType::LIST; break;
                        case 2: conf.BoxType = Detail::BoxGenerationType::VALIDATE; break;
                    }
                } 
                if(i == cur_mode) ImGui::SetItemDefaultFocus();
            }
            ImGui::EndCombo();
        }

        ImGui::Checkbox("Multithreading", &conf.MultiThreading);
        ImGui::Checkbox("Random Seed", &conf.UseRandomSeed);
        {
            seedval = conf.Seed;
            char buf[20];
            std::sprintf(buf, "%llu", seedval);
            if(ImGui::InputText("Seed", buf, IM_ARRAYSIZE(buf), ImGuiInputTextFlags_CallbackAlways, [](ImGuiInputTextCallbackData* data) {
                if (data->EventFlag == ImGuiInputTextFlags_CallbackAlways) {
                    uint64_t* p_value = (uint64_t*)data->UserData;
                    std::string input_text(data->Buf);
                    uint64_t value;
                    auto [ptr, ec] = std::from_chars(input_text.data(), input_text.data() + input_text.size(), value);
                    if (ec != std::errc::invalid_argument && ec != std::errc::result_out_of_range)
                        *p_value = value; 
                }
                return 0;
            }, &seedval)){
                if(!conf.UseRandomSeed) conf.Seed = seedval;
            }
        }
        ImGui::Checkbox("Allow Overlap", &conf.AllowOverlap);

        ImGui::Dummy({0, 4});
        ImGui::Separator();
        ImGui::Dummy({0, 4});

        ImGui::InputFloat2("Bounds", glm::value_ptr(conf.Bounds));
        ImGui::InputInt("Height", (int32_t*)&conf.Height);

        ImGui::Dummy({0, 4});
        ImGui::Separator();
        ImGui::Dummy({0, 4});

        ImGui::InputInt("Empty Tries", (int32_t*)&conf.EmptryTries);
        ImGui::InputInt("Max Empty Tries", (int32_t*)&conf.MaxEmptryTries);
        ImGui::InputInt("Look Ahead", (int32_t*)&conf.LookAheadSize);

        ImGui::Dummy({0, 4});
        ImGui::Separator();
        ImGui::Dummy({0, 4});

        if(conf.BoxType == Detail::BoxGenerationType::RANDOM){
            ImGui::Checkbox("Use Cubes", &conf.CubeRandomBoxes);
            ImGui::InputFloat("Min Box Volume", &conf.MinBoxVolume);
            ImGui::InputFloat("Max Box Volume", &conf.MaxBoxVolume);
        } else if(conf.BoxType == Detail::BoxGenerationType::LIST) {
            ImGui::Checkbox("Enforce Misses", &conf.EnforceMisses);
        } else if(conf.BoxType == Detail::BoxGenerationType::VALIDATE) {

            ImGui::InputInt("Box Count", (int32_t*)&conf.EvalBoxCount);
            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.f, 1.f, 0.f, 1.f));
            ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.0f, 0.0f, 0.0f, 1.0f));
            if(ImGui::Button("Evaluate\nProblem", ImVec2(100, 60))) {
                if(conf.UseRandomSeed){
                    conf.Seed = std::random_device()();
                }
                eval_res = Util::create_ground_truth(glm::vec<2, int32_t>{(int32_t)conf.Bounds.x, (int32_t)conf.Bounds.y}, conf.EvalBoxCount, conf.Seed);
            }

            ImGui::PopStyleColor(2);
            ImGui::SameLine();
            ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(1.f, 0.f, 0.0f, 1.f));
            if(ImGui::Button("Reset", ImVec2(100, 60))) {
                eval_res = nullptr;
            }
            ImGui::PopStyleColor();

            if(eval_res){
                const float ww = (ImGui::GetContentRegionAvail().x - 200) * 0.5f;
                const float ar = conf.Bounds.y / conf.Bounds.x;
                ImGui::Dummy({0, 10});
                ImGui::SetCursorPosX(ww);
                ImGui::Image((void*)&eval_res->tex_, 
                    ImVec2(200, int32_t(200.f * ar))
                );
            }
            
        }

        if(!ImGui::IsWindowHovered()) camera.update(GetFrameTime());
       
        ImGui::End();
        ImGui::Render();
        
        //------------------------

        BeginDrawing();
        ClearBackground(RAYWHITE);

        BeginMode3D(camera);
        DrawGrid(100, 5.0f);
        const auto ext = conf.Bounds * 0.5f;
        DrawCubeWiresV(Vector3{ext.y, conf.Height * 0.5f, ext.x}, Vector3{conf.Bounds.y, (float)conf.Height, conf.Bounds.x}, PURPLE);
        DrawLine3D(Vector3{0, 0, 0}, Vector3{(float)conf.Height, 0, 0}, RED);
        DrawLine3D(Vector3{0, 0, 0}, Vector3{0, (float)conf.Height, 0}, GREEN);
        DrawLine3D(Vector3{0, 0, 0}, Vector3{0, 0, (float)conf.Height}, BLUE);

        if(pr){
            if(cc++%30 == 0)
                b = pr->data_cpy();

            for (const auto& tr : b) {
                const Vector3 ns = {glm::length(glm::vec3(tr[0])), glm::length(glm::vec3(tr[1])), glm::length(glm::vec3(tr[2])) };
                const double temp = 1. - 1. / (conf.MaxBoxVolume - conf.MinBoxVolume) * ((ns.x * ns.y * ns.z) - conf.MinBoxVolume);
                const auto trans = glm::vec<3, float>(tr[3]);
                DrawCubeV(Vector3{ trans.x, trans.z, trans.y }, Vector3{ ns.x, ns.z, ns.y }, col(temp));
                DrawCubeWiresV(Vector3{ trans.x, trans.z, trans.y }, Vector3{ ns.x, ns.z, ns.y }, BLACK);
            }
        }

        EndMode3D();
      
        ImGui_ImplRaylib_RenderDrawData(ImGui::GetDrawData());
        DrawFPS(GetScreenWidth() - 100, 5);
        EndDrawing();

    }
    ImGui_ImplRaylib_Shutdown();
    CloseWindow();
    return 0;
}
