/****************************************************************************
** Meta object code from reading C++ file 'MeshSkeletonEditingDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../include/MeshSkeletonEditingDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshSkeletonEditingDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshSkeletonEditingDialog_t {
    QByteArrayData data[13];
    char stringdata0[159];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshSkeletonEditingDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshSkeletonEditingDialog_t qt_meta_stringdata_JMeshSkeletonEditingDialog = {
    {
QT_MOC_LITERAL(0, 0, 26), // "JMeshSkeletonEditingDialog"
QT_MOC_LITERAL(1, 27, 12), // "selectBranch"
QT_MOC_LITERAL(2, 40, 0), // ""
QT_MOC_LITERAL(3, 41, 13), // "displayBranch"
QT_MOC_LITERAL(4, 55, 11), // "getContours"
QT_MOC_LITERAL(5, 67, 13), // "deleteContour"
QT_MOC_LITERAL(6, 81, 12), // "removeBranch"
QT_MOC_LITERAL(7, 94, 11), // "splitBranch"
QT_MOC_LITERAL(8, 106, 14), // "reparameterize"
QT_MOC_LITERAL(9, 121, 9), // "enumNodes"
QT_MOC_LITERAL(10, 131, 6), // "getCap"
QT_MOC_LITERAL(11, 138, 8), // "getPoles"
QT_MOC_LITERAL(12, 147, 11) // "closeDialog"

    },
    "JMeshSkeletonEditingDialog\0selectBranch\0"
    "\0displayBranch\0getContours\0deleteContour\0"
    "removeBranch\0splitBranch\0reparameterize\0"
    "enumNodes\0getCap\0getPoles\0closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshSkeletonEditingDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      11,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   69,    2, 0x08 /* Private */,
       3,    0,   70,    2, 0x08 /* Private */,
       4,    0,   71,    2, 0x08 /* Private */,
       5,    0,   72,    2, 0x08 /* Private */,
       6,    0,   73,    2, 0x08 /* Private */,
       7,    0,   74,    2, 0x08 /* Private */,
       8,    0,   75,    2, 0x08 /* Private */,
       9,    0,   76,    2, 0x08 /* Private */,
      10,    0,   77,    2, 0x08 /* Private */,
      11,    0,   78,    2, 0x08 /* Private */,
      12,    0,   79,    2, 0x08 /* Private */,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JMeshSkeletonEditingDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshSkeletonEditingDialog *_t = static_cast<JMeshSkeletonEditingDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->selectBranch(); break;
        case 1: _t->displayBranch(); break;
        case 2: _t->getContours(); break;
        case 3: _t->deleteContour(); break;
        case 4: _t->removeBranch(); break;
        case 5: _t->splitBranch(); break;
        case 6: _t->reparameterize(); break;
        case 7: _t->enumNodes(); break;
        case 8: _t->getCap(); break;
        case 9: _t->getPoles(); break;
        case 10: _t->closeDialog(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshSkeletonEditingDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshSkeletonEditingDialog.data,
      qt_meta_data_JMeshSkeletonEditingDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshSkeletonEditingDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshSkeletonEditingDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshSkeletonEditingDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshSkeletonEditingDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshSkeletonEditingDialog"))
        return static_cast< Ui::MeshSkeletonEditingDialog*>(const_cast< JMeshSkeletonEditingDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshSkeletonEditingDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 11)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 11;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 11)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 11;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
