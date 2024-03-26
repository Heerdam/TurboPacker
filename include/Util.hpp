#pragma once

#include <format>
#include <string>

#include <raylib.h>
#include <raymath.h>

namespace Util {

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
