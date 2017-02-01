/****************************************************************************
** Meta object code from reading C++ file 'MeshRenderDialog.hpp'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.6.0)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "MeshRenderDialog.hpp"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'MeshRenderDialog.hpp' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.6.0. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_JMeshRenderDialog_t {
    QByteArrayData data[15];
    char stringdata0[219];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_JMeshRenderDialog_t, stringdata0) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_JMeshRenderDialog_t qt_meta_stringdata_JMeshRenderDialog = {
    {
QT_MOC_LITERAL(0, 0, 17), // "JMeshRenderDialog"
QT_MOC_LITERAL(1, 18, 18), // "openMaterialDialog"
QT_MOC_LITERAL(2, 37, 0), // ""
QT_MOC_LITERAL(3, 38, 18), // "openMeshlistDialog"
QT_MOC_LITERAL(4, 57, 18), // "openBoxColorDialog"
QT_MOC_LITERAL(5, 76, 12), // "setEnclosure"
QT_MOC_LITERAL(6, 89, 11), // "resetCamera"
QT_MOC_LITERAL(7, 101, 14), // "fitBoundSphere"
QT_MOC_LITERAL(8, 116, 11), // "fitBoundBox"
QT_MOC_LITERAL(9, 128, 15), // "openNodesDialog"
QT_MOC_LITERAL(10, 144, 15), // "openEdgesDialog"
QT_MOC_LITERAL(11, 160, 15), // "openFacesDialog"
QT_MOC_LITERAL(12, 176, 15), // "openCellsDialog"
QT_MOC_LITERAL(13, 192, 14), // "setPOVRayScene"
QT_MOC_LITERAL(14, 207, 11) // "closeDialog"

    },
    "JMeshRenderDialog\0openMaterialDialog\0"
    "\0openMeshlistDialog\0openBoxColorDialog\0"
    "setEnclosure\0resetCamera\0fitBoundSphere\0"
    "fitBoundBox\0openNodesDialog\0openEdgesDialog\0"
    "openFacesDialog\0openCellsDialog\0"
    "setPOVRayScene\0closeDialog"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_JMeshRenderDialog[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
      13,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: name, argc, parameters, tag, flags
       1,    0,   79,    2, 0x08 /* Private */,
       3,    0,   80,    2, 0x08 /* Private */,
       4,    0,   81,    2, 0x08 /* Private */,
       5,    0,   82,    2, 0x08 /* Private */,
       6,    0,   83,    2, 0x08 /* Private */,
       7,    0,   84,    2, 0x08 /* Private */,
       8,    0,   85,    2, 0x08 /* Private */,
       9,    0,   86,    2, 0x08 /* Private */,
      10,    0,   87,    2, 0x08 /* Private */,
      11,    0,   88,    2, 0x08 /* Private */,
      12,    0,   89,    2, 0x08 /* Private */,
      13,    0,   90,    2, 0x08 /* Private */,
      14,    0,   91,    2, 0x08 /* Private */,

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
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void JMeshRenderDialog::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        JMeshRenderDialog *_t = static_cast<JMeshRenderDialog *>(_o);
        Q_UNUSED(_t)
        switch (_id) {
        case 0: _t->openMaterialDialog(); break;
        case 1: _t->openMeshlistDialog(); break;
        case 2: _t->openBoxColorDialog(); break;
        case 3: _t->setEnclosure(); break;
        case 4: _t->resetCamera(); break;
        case 5: _t->fitBoundSphere(); break;
        case 6: _t->fitBoundBox(); break;
        case 7: _t->openNodesDialog(); break;
        case 8: _t->openEdgesDialog(); break;
        case 9: _t->openFacesDialog(); break;
        case 10: _t->openCellsDialog(); break;
        case 11: _t->setPOVRayScene(); break;
        case 12: _t->closeDialog(); break;
        default: ;
        }
    }
    Q_UNUSED(_a);
}

const QMetaObject JMeshRenderDialog::staticMetaObject = {
    { &QDialog::staticMetaObject, qt_meta_stringdata_JMeshRenderDialog.data,
      qt_meta_data_JMeshRenderDialog,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *JMeshRenderDialog::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *JMeshRenderDialog::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_JMeshRenderDialog.stringdata0))
        return static_cast<void*>(const_cast< JMeshRenderDialog*>(this));
    if (!strcmp(_clname, "Ui::MeshRenderDialog"))
        return static_cast< Ui::MeshRenderDialog*>(const_cast< JMeshRenderDialog*>(this));
    return QDialog::qt_metacast(_clname);
}

int JMeshRenderDialog::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QDialog::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 13)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 13;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 13)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 13;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
