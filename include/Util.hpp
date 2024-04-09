#pragma once

#include <format>
#include <string>
#include <filesystem>
#include <fstream>

#include <raylib.h>
#include <raymath.h>

#include <nlohmann/json.hpp>

#include <TP.hpp>

namespace Util {

    namespace Detail {

        struct EvalRes {
            glm::vec<2, int32_t> bounds_;
            std::vector<Color> colmap_;
            std::vector<std::pair<int32_t, TP::Detail::FBox<2, int32_t>>> decomp_;
            Texture2D tex_;
            TP::Detail::BoxList<float> boxlist_;
        };//EvalRes

    }//Detail

    [[nodiscard]] inline std::unique_ptr<Detail::EvalRes> create_ground_truth(
        const glm::vec<2, int32_t> _bounds,
        const int32_t _num_boxes,
        const uint64_t _seed
    ) {

        std::mt19937_64 rng (_seed);

        std::unique_ptr<Detail::EvalRes> out = std::make_unique<Detail::EvalRes>();
        out->bounds_ = _bounds;

        out->colmap_.resize(_num_boxes);
        for(int32_t i = 0; i < _num_boxes; ++i)
            out->colmap_[i] = ColorFromHSV(i * (360.f / (_num_boxes + 1)), 1.f, 1.f);

        std::vector<std::pair<int32_t, TP::Detail::Interval<float>>> bxs;
        bxs.resize(_num_boxes);

        const float median = out->bounds_.x * out->bounds_.y / _num_boxes;
        for(int32_t i = 0; i < _num_boxes; ++i)
            bxs[i] = { i, {median * 0.5f, median * 1.5f} };

        out->decomp_ = TP::squarify(TP::Detail::FBox<2, float>{ {0.f, 0.f}, { float(_bounds.x), float(_bounds.y) } }, bxs, rng());

        std::vector<unsigned char> img_data;
        img_data.resize(4 * _bounds.x * _bounds.y, 0);

        for(const auto[idx, b] : out->decomp_){

            //fill
            for(size_t y = b.min_.y; y < b.max_.y; ++y){
                for(size_t x = b.min_.x; x < b.max_.x; ++x){
                    const size_t i = x + y * _bounds.x;
                    if(4*i >= img_data.size()) continue;
                    img_data[4*i] = out->colmap_[idx].r;
                    img_data[4*i + 1] = out->colmap_[idx].g;
                    img_data[4*i + 2] = out->colmap_[idx].b;
                    img_data[4*i + 3] = 255;
                }
            }

            //border y
            for(size_t y = b.min_.y; y < b.max_.y; ++y){
                const size_t i1 = b.min_.x + y * _bounds.x;
                const size_t i2 = (b.max_.x - 1) + y * _bounds.x;
                if(4*i1 >= img_data.size()) continue;
                if(4*i2 >= img_data.size()) continue;
                img_data[4*i1] = 0;
                img_data[4*i1 + 1] = 0;
                img_data[4*i1 + 2] = 0;
                img_data[4*i1 + 3] = 255;

                img_data[4*i2] = 0;
                img_data[4*i2 + 1] = 0;
                img_data[4*i2 + 2] = 0;
                img_data[4*i2 + 3] = 255;
            }

            //border x
            for(size_t x = b.min_.x; x < b.max_.x; ++x){
                const size_t i1 = x + b.min_.y * _bounds.x;
                const size_t i2 = x + (b.max_.y - 1) * _bounds.x;
                if(4*i1 >= img_data.size()) continue;
                if(4*i2 >= img_data.size()) continue;
                img_data[4*i1] = 0;
                img_data[4*i1 + 1] = 0;
                img_data[4*i1 + 2] = 0;
                img_data[4*i1 + 3] = 255;

                img_data[4*i2] = 0;
                img_data[4*i2 + 1] = 0;
                img_data[4*i2 + 2] = 0;
                img_data[4*i2 + 3] = 255;
            }

        }

        Image img;
        img.mipmaps = 1;
        img.width = _bounds.x;
        img.height = _bounds.y;
        img.format = PixelFormat::PIXELFORMAT_UNCOMPRESSED_R8G8B8A8;
        img.data = img_data.data();

        out->tex_ = LoadTextureFromImage(img);

        return out;
    }//create_groundtruth

    template<class T = float>
    class Camera3d : public Camera3D {

        Vector2 cursorPos = GetMousePosition();
        bool isDirty = true;

    public:

        T camDist = 250.; 
        T rotAngle = 45; 
        T tiltAngle = 4; 
        T rotSpeed = 0.5; 
        T moveSpeed = 50.; 
        T zoomSpeed = 15.;

        Camera3d() {
            projection = CAMERA_PERSPECTIVE;
            fovy = 65.f;
        }

        void update(const T _dt) {

            if (IsMouseButtonDown(1)) {
                Vector2 newPos = GetMousePosition();

                rotAngle -= (newPos.x - cursorPos.x) * rotSpeed;
                tiltAngle -= (newPos.y - cursorPos.y) * rotSpeed;

                if (tiltAngle > 89)
                    tiltAngle = 89;
                if (tiltAngle < -89)
                    tiltAngle = -89;

                isDirty = true;
            }

            cursorPos = GetMousePosition();

            Vector3 moveVec = { 0, 0, 0 };

            if (IsKeyDown(KEY_W))
                moveVec.z = -moveSpeed * _dt;
            if (IsKeyDown(KEY_S))
                moveVec.z = moveSpeed * _dt;

            if (IsKeyDown(KEY_A))
                moveVec.x = -moveSpeed * _dt;
            if (IsKeyDown(KEY_D))
                moveVec.x = moveSpeed * _dt;
        
            camDist -= zoomSpeed * GetMouseWheelMove();
            if (camDist < 1)
                camDist = 1;

            Vector3 camPos = { 0, 0, camDist};
            Matrix tiltMat = MatrixRotateX(tiltAngle * DEG2RAD);
            Matrix rotMat = MatrixRotateY(rotAngle * DEG2RAD); 
            Matrix mat = MatrixMultiply(tiltMat, rotMat); 

            camPos = Vector3Transform(camPos, mat);   
            moveVec = Vector3Transform(moveVec, rotMat); 
            target = Vector3Add(target, moveVec); 
            position = Vector3Add(target, camPos); 

            isDirty = false;
        }


    };

}//Util

namespace Disk {

    using json = nlohmann::json;

    template<class T, template<typename> class COSTFUNCTION, class HEIGHTMAP_T, class R_T, uint32_t BUCKET_SIZE, class HEIGHTMAP_ALLOCATOR>
    void load(const std::filesystem::path& _path) {

        using Conf = TP::Config<T, COSTFUNCTION, HEIGHTMAP_T, R_T, BUCKET_SIZE, HEIGHTMAP_ALLOCATOR>;

        std::ifstream i(_path);
        if(!i.good()) return Conf();

        Conf o;



        return o;

    }//load

    template<class T, template<typename> class CF, class H_T, class R_T, uint32_t BS, class HA>
    void save(
        const TP::Config<T, CF, H_T, R_T, BS, HA>& _conf,
        const std::filesystem::path& _path
    ) {

        json o;
        o["MultiThreading"] = _conf.MultiThreading;
        o["NumThreads"] = _conf.NumThreads;
        o["UseRandomSeed"] = _conf.UseRandomSeed;
        o["Seed"] = _conf.Seed;
        o["Bounds"] = { {"x", _conf.Bounds.x }, {"y", _conf.Bounds.y } };
        o["Height"] = _conf.Height;
        o["BoxType"] = int32_t(_conf.BoxType);
        o["AllowOverlap"] = _conf.AllowOverlap;
        o["EmptryTries"] = _conf.EmptryTries;
        o["LookAheadSize"] = _conf.LookAheadSize;
        o["AllowedPermutations"] = _conf.AllowedPermutations;
        o["CubeRandomBoxes"] = _conf.CubeRandomBoxes;
        o["MinBoxVolume"] = _conf.MinBoxVolume;
        o["MaxBoxVolume"] = _conf.MaxBoxVolume;
        json bl = json::object();
        bl["list_"] = json::array();
        bl["shuffle_boxes_"] = l.shuffle_boxes_;
        for(const TP::Detail::BoxEntry<T>& l : _conf.BoxList.list_){
            json be = json::object();
            be["size_"] = { { "x", l.size_.x }, { "y", l.size_.y }, { "z", l.size_.z } };
            be["count_"] = l.count_;
            be["perms_"] = l.perms_;
            bl["list_"].push_back(be);    
        }
        o["BoxList"] = bl;
        o["EnforceMisses"] = _conf.EnforceMisses;
        o["EvalBoxCount"] = _conf.EvalBoxCount;

        std::ofstream s(_path);
        s << std::setw(4) << o;
    }//save

}//Json
